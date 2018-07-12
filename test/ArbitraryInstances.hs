{-# LANGUAGE FlexibleInstances, RankNTypes, KindSignatures, DataKinds #-}
{-# OPTIONS_GHC -fno-warn-orphans #-}

-- | Orphan "Arbitrary" and related instances for testing purposes. 

module ArbitraryInstances where

import Control.Applicative
import Control.Lens((^.))
import Control.Monad
import Geodetics.Altitude
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.GridScale
import Geodetics.Ellipsoids
import Geodetics.Latitude
import Geodetics.Longitude
import Geodetics.Path
import Geodetics.Stereographic as SG
import Geodetics.TransverseMercator as TM
import Numeric.Units.Dimensional.Prelude
import qualified Prelude ()
import Test.QuickCheck



-- | Shrink using a dimension, so that shrunk values are round numbers in that dimension.
shrinkDimension :: forall a (d :: Dimension) (m :: Metricality) .
                   (Fractional a, Arbitrary a) => Unit m d a -> Quantity d a -> [Quantity d a]
shrinkDimension u v = (*~ u) <$> shrink (v /~ u)

-- | Wrapper for arbitrary angles.
newtype Bearing = Bearing (Dimensionless Double)

instance Show Bearing where
   show (Bearing b) = "Bearing " ++ showAngle b

instance Arbitrary Bearing where
   arbitrary = Bearing <$> (*~ degree) <$> choose (-180,180)
   shrink (Bearing b) = Bearing <$> shrinkDimension degree b
   
   
newtype Azimuth = Azimuth (Dimensionless Double)

instance Show Azimuth where
   show (Azimuth a) = "Azimuth " ++ showAngle a
   
instance Arbitrary Azimuth where
   arbitrary = Azimuth <$> (*~ degree) <$> choose (0,90)
   shrink (Azimuth a) = Azimuth <$> shrinkDimension degree a
   
   
-- | Wrapper for arbitrary distances up to 10,000 km
newtype Distance = Distance (Length Double) deriving (Show)

instance Arbitrary Distance where
   arbitrary = Distance <$> (*~ kilo meter) <$> choose (0,10000)
   shrink (Distance d) = Distance <$> shrinkDimension (kilo meter) d
   

-- | Wrapper for arbitrary distances up to 1,000 km
newtype Distance2 = Distance2 (Length Double) deriving (Show)

instance Arbitrary Distance2 where
   arbitrary = Distance2 <$> (*~ kilo meter) <$> choose (0,1000)
   shrink (Distance2 d) = Distance2 <$> shrinkDimension (kilo meter) d

-- | Wrapper for arbitrary altitudes up to 10 km
newtype Altitude = Altitude (Length Double) deriving (Show)

instance Arbitrary Altitude where
   arbitrary = Altitude <$> (*~ kilo meter) <$> choose (0,10)
   shrink (Altitude h) = Altitude <$> shrinkDimension (kilo meter) h


-- | Wrapper for arbitrary dimensionless numbers (-10 .. 10)
newtype Scalar = Scalar (Dimensionless Double) deriving (Show)

instance Arbitrary Scalar where
   arbitrary = Scalar <$> (*~ one) <$> choose (-10,10)
   shrink (Scalar s) = Scalar <$> shrinkDimension one s


-- | Wrapper for arbitrary grid references.
newtype GridRef = GridRef String deriving Show

instance Arbitrary GridRef where
   arbitrary = do
      n <- choose (0,4)
      c1 <- elements "HJNOST" -- General vicinity of UK
      c2 <- elements $ ['A'..'H'] ++ ['J'..'Z']
      dx <- vectorOf n $ choose ('0','9')
      dy <- vectorOf n $ choose ('0','9')
      return $ GridRef $ c1 : c2 : (dx ++ dy)
   shrink = shrinkNothing


-- | Generate in range +/- <arg> m.
genOffset :: Double -> Gen (Length Double)
genOffset d = (*~ meter) <$> choose (-d, d)

genAlt :: Gen (Length Double)
genAlt = (*~ meter) <$> choose (0,10000)


genLatitude :: Gen (Dimensionless Double)
genLatitude = (*~ degree) <$> choose (-90,90)

genLongitude :: Gen (Dimensionless Double)
genLongitude = (*~ degree) <$> choose (-180,180)

genSeconds :: Gen (Dimensionless Double)
genSeconds = (*~ arcsecond) <$> choose (-10,10)
    

-- | Shrinking with the original value preserved. Used for shrinking records.  See 
-- http://stackoverflow.com/questions/14006005/idiomatic-way-to-shrink-a-record-in-quickcheck for details.
shrink' :: (Arbitrary a) => a -> [a]
shrink' x = x : shrink x

-- | Shrink a quantity in the given units.
shrinkQuantity :: forall a (d :: Dimension) (m :: Metricality).
                  (Arbitrary a, Fractional a) => Unit m d a -> Quantity d a -> [Quantity d a]
shrinkQuantity u q = map (*~ u) $ shrink' $ q /~ u

shrinkLength :: (Arbitrary a, Fractional a) => Length a -> [Length a]
shrinkLength = shrinkQuantity meter

shrinkUnit :: (Arbitrary a, Fractional a) => Dimensionless a -> [Dimensionless a]
shrinkUnit = shrinkQuantity one

shrinkAngle :: (Arbitrary a, Floating a) => Dimensionless a -> [Dimensionless a]
shrinkAngle = shrinkQuantity degree


instance Arbitrary Helmert where
   arbitrary = 
      Helmert <$> genOffset 300 <*> genOffset 300 <*> genOffset 300 <*> 
         ((*~ one) <$> choose (-5,10)) <*>
         genSeconds <*> genSeconds <*> genSeconds
   shrink h = 
      tail $ Helmert <$> shrinkLength (h ^. cX) <*> shrinkLength (h ^. cY) <*> shrinkLength (h ^. cZ) <*>
         shrinkUnit (h ^. helmertScale) <*>
         shrinkUnit (h ^. rX) <*> shrinkUnit (h ^. rY) <*> shrinkUnit (h ^. rZ)

instance Arbitrary TRF where
   arbitrary =
      TRF <$>
         ((*~ meter) <$> choose (6378100, 6378400)) <*>                  -- majorRadius
         ((*~ one) <$> choose (297,300)) <*>                             -- flatR
         arbitrary                                                       -- helmert
   shrink e = tail $ TRF (e ^. majorRadius) (e ^. flatR) <$> shrink' (e ^. helmert)
        
instance Arbitrary Geodetic where
   arbitrary = 
      Geodetic <$> genLatitude <*> genLongitude <*> genOffset 1 <*> arbitrary
   shrink g = 
      tail $
         Geodetic <$>
         shrinkAngle (g ^. latitude) <*>
         shrinkAngle (g ^. longitude) <*> 
         shrinkLength (g ^. altitude) <*>
         shrink' (g ^. trf)

instance Arbitrary (GridPoint GridTM) where
   arbitrary = GridPoint <$> genOffset 100000 <*> genOffset 100000 <*> genOffset 1 <*> arbitrary
   shrink p = 
      tail $ GridPoint <$> 
         shrinkLength (p ^. eastings) <*> 
         shrinkLength (p ^. northings) <*> 
         shrinkLength (p ^. altitude) <*> 
         shrink' (p ^. gridBasis)


instance Arbitrary (GridPoint GridStereo) where
   arbitrary = GridPoint <$> genOffset 100000 <*> genOffset 100000 <*> genOffset 1 <*> arbitrary
   shrink p =
      tail $ GridPoint <$> 
         shrinkLength (p ^. eastings) <*> 
         shrinkLength (p ^. northings) <*> 
         shrinkLength (p ^. altitude) <*> 
         shrink' (p ^. gridBasis)


instance Arbitrary GridTM where
   arbitrary = mkGridTM <$> arbitrary <*> arbitrary <*> ((*~ one) <$> choose (0.95,1.0))
   shrink tm = tail $ mkGridTM <$> shrink' (tm ^. geodetic) <*> shrink' (tm ^. gridOffsetL) <*> [tm ^. gridScale]
   
   
instance Arbitrary GridOffset where
   arbitrary = GridOffset <$> genOffset 100000 <*> genOffset 100000 <*> genAlt
   shrink d = tail $ GridOffset <$> shrinkLength (d ^. deltaEastL) <*> shrinkLength (d ^. deltaNorthL) <*> shrinkLength (d ^. deltaAltitudeL)


instance Arbitrary GridStereo where
   arbitrary = mkGridStereo <$> arbitrary <*> arbitrary <*> ((*~ one) <$> choose (0.95,1.0))
   shrink sg = tail $ mkGridStereo <$> shrink' (sg ^. geodetic) <*> shrink' (sg ^. gridOffsetL) <*> [sg ^. gridScale]
   

-- | Wrapper for arbitrary rays, along with creation parameters for printing and shrinking.
data Ray = Ray Geodetic (Angle Double) (Angle Double)

instance Show Ray where
   show (Ray p0 b e ) = "(Ray " ++ show p0 ++ ", " ++ showAngle b ++ ", " ++ showAngle e ++ ")"

getRay :: Ray -> Path
getRay (Ray p0 b e) = rayPath p0 b e

instance Arbitrary Ray where
   arbitrary = do
      p0 <- arbitrary
      b <- (*~ degree) <$> choose (-180,180)
      e <- (*~ degree) <$> choose (0,90)
      return $ Ray p0 b e
   shrink (Ray p0 b e) = tail $ do
      p0' <- shrink' p0
      b' <- shrinkAngle b
      e' <- shrinkAngle e
      return $ Ray p0' b' e'
      
     
-- | Two rhumb paths starting not more than 1000 km apart.
data RhumbPaths2 = RP2 {
      rp2Point0 :: Geodetic,
      rp2Bearing0 :: Bearing,
      rp2Distance :: Distance2,
      rp2Bearing1 :: Bearing,
      rp2Bearing2 :: Bearing
   }
   
instance Show RhumbPaths2 where
   show rp2 = show (pt1, Bearing b1) ++ show (pt2, Bearing b2)
      where 
         (p1, p2) = mk2RhumbPaths rp2
         (pt1, b1, _) = pathFunc p1 (0 *~ meter)
         (pt2, b2, _) = pathFunc p2 (0 *~ meter)
          
instance Arbitrary RhumbPaths2 where
   arbitrary = RP2 
      <$> arbitrary `suchThat` ((< 70 *~ degree) . abs . (^. latitude))
      <*> arbitrary <*> arbitrary <*> arbitrary <*> arbitrary
   shrink rp = 
      tail $ RP2 <$> 
         shrink' (rp2Point0 rp) <*> 
         shrink' (rp2Bearing0 rp) <*> 
         shrink' (rp2Distance rp) <*> 
         shrink' (rp2Bearing1 rp) <*> 
         shrink' (rp2Bearing2 rp)

mk2RhumbPaths :: RhumbPaths2 -> (Path, Path)
mk2RhumbPaths (RP2 pt0 (Bearing b0) (Distance2 d) (Bearing b1) (Bearing b2)) =
   (path1, path2)
   where
      path0 = rhumbPath pt0 b0
      path1 = rhumbPath pt0 b1
      (pt2, _, _) = pathFunc path0 d
      path2 = rhumbPath pt2 b2
 
