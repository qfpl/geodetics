{-# LANGUAGE NoImplicitPrelude #-}

{- |
The following is based on equations in Section 1.4.7.1 in 
OGP Surveying and Positioning Guidance Note number 7, part 2 – August 2006
<http://ftp.stu.edu.tw/BSD/NetBSD/pkgsrc/distfiles/epsg-6.11/G7-2.pdf>
-}

{-# LANGUAGE MultiParamTypeClasses, FlexibleInstances #-}
module Geodetics.Stereographic (
   GridStereo(..),
   HasGridStereo(..),
   mkGridStereo
) where

import Control.Lens(Lens', (^.))
import Geodetics.Altitude(HasAltitude(altitude))
import Geodetics.Ellipsoids
import Geodetics.Geodetic(Geodetic(Geodetic), HasGeodetic(geodetic))
import Geodetics.Grid(HasGridOffset(gridOffsetL), GridClass(gridEllipsoid), GridOffset, GridPoint(GridPoint), toGrid, fromGrid, applyOffset, offsetNegate, gridBasis)
import Geodetics.GridScale(HasGridScale(gridScale))
import Geodetics.Latitude(HasLatitude(latitudeL))
import Geodetics.Longitude(HasLongitude(longitudeL))
import Numeric.Units.Dimensional.Prelude


-- | A stereographic projection with its origin at an arbitrary point on Earth, other than the poles.
data GridStereo = GridStereo
      Geodetic -- ^ Point where the plane of projection touches the ellipsoid. Often known as the Natural Origin.
      GridOffset  -- ^ Grid position of the tangent point. Often known as the False Origin.
      (Dimensionless Double) -- ^ Scaling factor that balances the distortion between the center and the edges. 
                                         -- Should be slightly less than unity.
      
      -- Memoised parameters derived from the tangent point.
      (Length Double)
      (Dimensionless Double)
      (Dimensionless Double)
      (Dimensionless Double)
      (Dimensionless Double)
      (Angle Double)
      (Length Double)
      (Length Double)
   deriving (Show)

class HasGridStereo a where
   gridStereo ::
     Lens' a GridStereo
   gridRL ::
     Lens' a (Length Double)
   gridNL :: 
     Lens' a (Dimensionless Double)
   gridCL :: 
     Lens' a (Dimensionless Double)
   gridSinL :: 
     Lens' a (Dimensionless Double)
   gridCosL :: 
     Lens' a (Dimensionless Double)
   gridLatCL :: 
     Lens' a (Angle Double)
   gridGL :: 
     Lens' a (Length Double)
   gridHL :: 
     Lens' a (Length Double)
   gridRL =
      gridStereo . gridRL
   {-# INLINE gridRL #-}
   gridNL =
      gridStereo . gridNL
   {-# INLINE gridNL #-}
   gridCL =
      gridStereo . gridCL
   {-# INLINE gridCL #-}
   gridSinL =
      gridStereo . gridSinL
   {-# INLINE gridSinL #-}
   gridCosL =
      gridStereo . gridCosL
   {-# INLINE gridCosL #-}
   gridLatCL =
      gridStereo . gridLatCL
   {-# INLINE gridLatCL #-}
   gridGL =
      gridStereo . gridGL
   {-# INLINE gridGL #-}
   gridHL =
      gridStereo . gridHL
   {-# INLINE gridHL #-}

instance HasGridScale GridStereo where
   gridScale k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o x r dn dc dsin dcos c g h) (k s)
   {-# INLINE gridScale #-}

instance HasGridStereo GridStereo where
   gridStereo =
      id
   gridRL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s x dn dc dsin dcos c g h) (k r)
   {-# INLINE gridRL #-}
   gridNL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r x dc dsin dcos c g h) (k dn)
   {-# INLINE gridNL #-}
   gridCL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r dn x dsin dcos c g h) (k dc)
   {-# INLINE gridCL #-}
   gridSinL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r dn dc x dcos c g h) (k dsin)
   {-# INLINE gridSinL #-}
   gridCosL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r dn dc dsin x c g h) (k dcos)
   {-# INLINE gridCosL #-}
   gridLatCL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r dn dc dsin dcos x g h) (k c)
   {-# INLINE gridLatCL #-}
   gridGL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r dn dc dsin dcos c x h) (k g)
   {-# INLINE gridGL #-}
   gridHL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t o s r dn dc dsin dcos c g x) (k h)
   {-# INLINE gridHL #-}

instance HasGeodetic GridStereo where
   geodetic k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo x o s r dn dc dsin dcos c g h) (k t)

instance HasAltitude GridStereo where
   altitude =
      geodetic . altitude

instance HasLatitude GridStereo where
   latitudeL =
      geodetic . latitudeL

instance HasLongitude GridStereo where
   longitudeL =
      geodetic . longitudeL

instance HasEllipsoid GridStereo where
   ellipsoid =
      geodetic . ellipsoid

instance HasGridOffset GridStereo where
   gridOffsetL k (GridStereo t o s r dn dc dsin dcos c g h) =
      fmap (\x -> GridStereo t x s r dn dc dsin dcos c g h) (k o)

-- | Create a stereographic projection. The tangency point must not be one of the poles.  
mkGridStereo :: Geodetic -> GridOffset -> Dimensionless Double -> GridStereo
mkGridStereo tangent origin scale = GridStereo tangent origin scale r n c sinLatC1 (sqrt $ _1 - sinLatC1 * sinLatC1) (asin sinLatC1) g h
   where 
      -- The reference seems to use χO to refer to two slightly different values. 
      -- Here these will be called LatC0 and LatC1.
      ellipse = tangent ^. ellipsoid
      op :: Num a => Quantity d a -> Quantity d a    -- Values of longitude, tangent longitude, E and N
      op = if tangent ^. latitudeL < _0 then negate else id  -- must be negated in the southern hemisphere.
      lat0 = op $ (tangent ^. latitudeL)
      sinLat0 = sin lat0
      e2 = eccentricity2 ellipse
      e = sqrt e2
      r = sqrt $ meridianRadius ellipse lat0 * primeVerticalRadius ellipse lat0
      n = sqrt $ _1 + ((e2 * cos lat0 ^ pos4)/(_1 - e2))
      s1 = (_1 + sinLat0) / (_1 - sinLat0)
      s2 = (_1 - e * sinLat0) / (_1 + e * sinLat0)
      w1 = (s1 * s2 ** e) ** n
      sinLatC0 = (w1 - _1)/(w1 + _1)
      c = ((n + sin lat0) * (_1 - sinLatC0)) / ((n - sin lat0) * (_1 + sinLatC0))
      w2 = c * w1
      sinLatC1 = (w2 - _1)/(w2 + _1)
      g = _2 * r * scale * tan (pi/_4 - latC1/_2)
      h = _4 * r * scale * tan latC1 + g
      latC1 = asin sinLatC1
      

instance GridClass GridStereo where
   toGrid grid geo = applyOffset (grid ^. gridOffsetL) $ GridPoint east north (geo ^. altitude) grid
      where
         op :: Num a => Quantity d a -> Quantity d a    -- Values of longitude, tangent longitude, E and N
         op = if (grid ^. latitudeL) < _0 then negate else id  -- must be negated in the southern hemisphere.
         sinLatC = (w - _1)/(w + _1)
         cosLatC = sqrt $ _1 - sinLatC * sinLatC
         longC = (grid ^. gridNL) * (op (geo ^. longitudeL) - long0) + long0
         w = (grid ^. gridCL) * (sA * sB ** e) ** (grid ^. gridNL)
         sA = (_1+sinLat) / (_1 - sinLat)
         sB = (_1 - e*sinLat) / (_1 + e*sinLat)
         sinLat = sin $ op $ (geo ^. latitudeL)
         e = sqrt $ eccentricity2 $ (geo ^. ellipsoid)
         long0 = op $ (grid ^. longitudeL)
         b = _1 + sinLatC * (grid ^. gridSinL) + cosLatC * (grid ^. gridCosL) * cos (longC - long0)
         east = _2 * (grid ^. gridRL) * (grid ^. gridScale) * cosLatC * sin (longC - long0) / b
         north = _2 * (grid ^. gridRL) * (grid ^. gridScale) * (sinLatC * (grid ^. gridCosL) - cosLatC * (grid ^. gridSinL) * cos (longC - long0)) / b
   
   fromGrid gp = 
      {- trace (    -- Remove comment brackets for debugging.
         "fromGrid values:\n   i = " ++ show i ++ "\n   j = " ++ show j ++
         "\n   longC = " ++ show longC ++ "\n   long = " ++ show long ++
         "\n   latC = " ++ show latC ++
         "\n   lat1 = " ++ show lat1 ++ "\n   latN = " ++ show latN ) $ -}
         Geodetic (op latN) (op long) height $ gridEllipsoid grid
      where
         op :: Num a => Quantity d a -> Quantity d a                   -- Values of longitude, tangent longitude, E and N
         op = if (grid ^. latitudeL) < _0 then negate else id  -- must be negated in the southern hemisphere.
         GridPoint east north height _ = applyOffset (offsetNegate $ grid ^. gridOffsetL) gp
         east' = east
         north' = north
         grid = gp ^. gridBasis
         long0 = op (grid ^. longitudeL)
         i = atan2 east' ((grid ^. gridHL) + north')
         j = atan2 east' ((grid ^. gridGL) - north') - i
         latC = (grid ^. gridLatCL) + _2 * atan2 (north' - east' * tan (j/_2)) (_2 * (grid ^. gridRL) * (grid ^. gridScale))
         longC = j + _2 * i + long0
         sinLatC = sin latC
         long = (longC - long0) / (grid ^. gridNL) + long0
         isoLat = log ((_1 + sinLatC) / ((grid ^. gridCL) * (_1 - sinLatC))) / (_2 * (grid ^. gridNL))
         lat1 = _2 * atan (exp isoLat) - pi/_2
         next lat = lat - (isoN - isoLat) * cos lat * (_1 - e2 * sin lat ^ pos2) / (_1 - e2)
            where isoN = isometricLatitude (gridEllipsoid grid) lat
                  e2 = eccentricity2 $ gridEllipsoid grid
         lats = iterate next lat1
         latN = snd . head . dropWhile (\(v1, v2) -> abs (v1-v2) > 0.01 *~ arcsecond) . zip lats . tail $ lats 
            
   gridEllipsoid = (^. ellipsoid)
