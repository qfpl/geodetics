{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE MultiParamTypeClasses, FunctionalDependencies, FlexibleInstances #-}

module Geodetics.Grid (
   -- ** Grid types
   GridClass (..),
   GridPoint (..),
   HasGridPoint(..),
   GridOffset (..),
   HasGridOffset (..),
   -- ** Grid operations
   polarOffset,
   offsetScale,
   offsetNegate,
   applyOffset,
   offsetDistance,
   offsetDistanceSq,
   offsetBearing,
   gridOffset,
   -- ** Unsafe conversion
   unsafeGridCoerce,
   -- ** Utility functions for grid references
   fromGridDigits,
   toGridDigits
) where

import Control.Lens(Lens', (^.), _Wrapped')
import Data.Char
import Data.Function
import Data.Monoid (Monoid)
import Data.Semigroup (Semigroup, (<>))
import Geodetics.Geodetic
import Geodetics.Types.Altitude
import Geodetics.Types.TRF
import Numeric.Units.Dimensional.Prelude hiding ((.))
import qualified Prelude as P

-- | A Grid is a two-dimensional projection of the ellipsoid onto a plane. Any given type of grid can
-- usually be instantiated with parameters such as a tangential point or line, and these parameters
-- will include the terrestrial reference frame ("Ellipsoid" in this library) used as a foundation. 
-- Hence conversion from a geodetic to a grid point requires the \"basis\" for the grid in question,
-- and grid points carry that basis with them because without it there is no defined relationship
-- between the grid points and terrestrial positions.
class GridClass r where
   fromGrid :: GridPoint r -> Geodetic
   toGrid :: r -> Geodetic -> GridPoint r
   gridTRF :: r -> TRF

-- | A point on the specified grid. 
data GridPoint r = GridPoint 
   (Length Double)
   (Length Double)
   (Length Double)
   r
   deriving (Show)

class HasGridPoint a r | a -> r where
   gridPoint ::
      Lens' a (GridPoint r)
   eastings ::
      Lens' a (Length Double)
   {-# INLINE eastings #-}
   northings ::
      Lens' a (Length Double)
   {-# INLINE northings #-}
   altGP ::
      Lens' a (Length Double)
   {-# INLINE altGP #-}
   gridBasis ::
      Lens' a r
   {-# INLINE gridBasis #-}
   altGP =
      gridPoint . altGP
   eastings =
      gridPoint . eastings
   gridBasis =
      gridPoint . gridBasis
   northings =
      gridPoint . northings

instance HasGridPoint (GridPoint r) r where
   {-# INLINE eastings #-}
   {-# INLINE northings #-}
   {-# INLINE altGP #-}
   {-# INLINE gridBasis #-}
   gridPoint = P.id
   eastings k (GridPoint e n a g) =
      fmap (\x -> GridPoint x n a g) (k e)
   northings k (GridPoint e n a g) =
      fmap (\x -> GridPoint e x a g) (k n)
   altGP k (GridPoint e n a g) =
      fmap (\x -> GridPoint e n x g) (k a)
   gridBasis k (GridPoint e n a g) =
      fmap (\x -> GridPoint e n a x) (k g)

instance Eq (GridPoint r) where
   p1 == p2  = 
      (p1 ^. eastings) == (p2 ^. eastings) && 
      (p1 ^. northings) == (p2 ^. northings) && 
      (p1 ^. altGP) == (p2 ^. altGP)

instance HasAltitude (GridPoint g) where
   altitude k (GridPoint e n a b) = 
      fmap (\(Altitude x) -> GridPoint e n x b) (k (Altitude a))



-- | A vector relative to a point on a grid.
-- Operations that use offsets will only give
-- meaningful results if all the points come from the same grid.
-- 
-- The monoid instance is the sum of offsets.
data GridOffset = GridOffset
  (Length Double)
  (Length Double)
  (Length Double)
  deriving (Eq, Show)

class HasGridOffset a where
   gridOffsetL ::
     Lens' a GridOffset
   deltaEastL ::
     Lens' a (Length Double)
   {-# INLINE deltaEastL #-}
   deltaNorthL ::
     Lens' a (Length Double)
   {-# INLINE deltaNorthL #-}
   deltaAltitudeL ::
     Lens' a (Length Double)
   {-# INLINE deltaAltitudeL #-}
   deltaEastL =
      gridOffsetL . deltaEastL
   deltaNorthL =
      gridOffsetL . deltaNorthL
   deltaAltitudeL =
      gridOffsetL . deltaAltitudeL

instance HasGridOffset GridOffset where
   {-# INLINE deltaEastL #-}
   {-# INLINE deltaNorthL #-}
   {-# INLINE deltaAltitudeL #-}
   gridOffsetL =
      P.id
   deltaEastL k (GridOffset e n a) =
      fmap (\x -> GridOffset x n a) (k e)
   deltaNorthL k (GridOffset e n a) =
      fmap (\x -> GridOffset e x a) (k n)
   deltaAltitudeL k (GridOffset e n a) =
      fmap (\x -> GridOffset e n x) (k a)

instance Semigroup GridOffset where
  g1 <> g2 = GridOffset (g1 ^. deltaEastL + g2 ^. deltaEastL)
                        (g1 ^. deltaNorthL + g2 ^. deltaNorthL)
                        (g1 ^. deltaAltitudeL + g2 ^. deltaAltitudeL)

instance Monoid GridOffset where
   mempty = GridOffset _0 _0 _0
   mappend = (<>)

-- | An offset defined by a distance and a bearing to the right of North.
--
-- There is no elevation parameter because we are using a plane to approximate an ellipsoid,
-- so elevation would not provide a useful result.  If you want to work with elevations
-- then "rayPath" will give meaningful results.
polarOffset :: Length Double -> Angle Double -> GridOffset
polarOffset r d = GridOffset (r * sin d) (r * cos d) _0


-- | Scale an offset by a scalar.
offsetScale :: Dimensionless Double -> GridOffset -> GridOffset
offsetScale s off = GridOffset (off ^. deltaEastL * s)
                               (off ^. deltaNorthL * s)
                               (off ^. deltaAltitudeL * s)

-- | Invert an offset.
offsetNegate :: GridOffset -> GridOffset
offsetNegate off = GridOffset (negate $ off ^. deltaEastL)
                              (negate $ off ^. deltaNorthL)
                              (negate $ off ^. deltaAltitudeL)


-- Add an offset on to a point to get another point.
applyOffset :: GridOffset -> GridPoint g -> GridPoint g
applyOffset off p = GridPoint ((p ^. eastings) + off ^. deltaEastL)
                           ((p ^. northings) + off ^. deltaNorthL)
                           ((p ^. altitude . _Wrapped') + off ^. deltaAltitudeL)
                           (p ^. gridBasis)


-- | The distance represented by an offset.
offsetDistance :: GridOffset -> Length Double
offsetDistance = sqrt . offsetDistanceSq


-- | The square of the distance represented by an offset.
offsetDistanceSq :: GridOffset -> Area Double
offsetDistanceSq off = 
   (off ^. deltaEastL) ^ pos2 + (off ^. deltaNorthL) ^ pos2 + (off ^. deltaAltitudeL) ^ pos2

              
-- | The direction represented by an offset, as bearing to the right of North.
offsetBearing :: GridOffset -> Angle Double
offsetBearing off = atan2 (off ^. deltaEastL) (off ^. deltaNorthL)


-- | The offset required to move from p1 to p2.             
gridOffset :: GridPoint g -> GridPoint g -> GridOffset
gridOffset p1 p2 = GridOffset ((p2 ^. eastings) - (p1 ^. eastings))
                              ((p2 ^. northings) - (p1 ^. northings))
                              ((p2 ^. altitude . _Wrapped') - (p1 ^. altitude . _Wrapped'))


-- | Coerce a grid point of one type into a grid point of a different type, 
-- but with the same easting, northing and altitude. This is unsafe because it
-- will produce a different position unless the two grids are actually equal.
--
-- It should be used only to convert between distinguished grids (e.g. "UkNationalGrid") and
-- their equivalent numerical definitions.
unsafeGridCoerce :: b -> GridPoint a -> GridPoint b
unsafeGridCoerce base p = GridPoint (p ^. eastings) (p ^. northings) (p ^. altitude . _Wrapped') base



-- | Convert a list of digits to a distance. The first argument is the size of the
-- grid square within which these digits specify a position. The first digit is
-- in units of one tenth of the grid square, the second one hundredth, and so on.
-- The first result is the lower limit of the result, and the second is the size
-- of the specified offset.
-- So for instance @fromGridDigits (100 *~ kilo meter) "237"@ will return
--
-- > Just (23700 meters, 100 meters)
--
-- If there are any non-digits in the string then the function returns @Nothing@.
fromGridDigits :: Length Double -> String -> Maybe (Length Double, Length Double)
fromGridDigits sq ds = if all isDigit ds then Just (d, p) else Nothing
   where
      n = length ds
      d = sum $ zipWith (*) 
         (map ((*~ one) . fromIntegral . digitToInt) ds) 
         (tail $ iterate (/ (10 *~ one)) sq)
      p = sq / ((10 *~ one) ** (fromIntegral n *~ one))
      
-- | Convert a distance into a digit string suitable for printing as part
-- of a grid reference. The result is the nearest position to the specified
-- number of digits, expressed as an integer count of squares and a string of digits.
-- If any arguments are invalid then @Nothing@ is returned.
toGridDigits ::
   Length Double    -- ^ Size of enclosing grid square. Must be at least 1 km.
   -> Int           -- ^ Number of digits to return. Must be positive.
   -> Length Double -- ^ Offset to convert into grid.
   -> Maybe (Integer, String)
toGridDigits sq n d =
   if sq < (1 *~ kilo meter) || n < 0 || d < _0 
   then Nothing
   else
      Just (sqs, pad)
   where
      p :: Integer
      p = 10 P.^ n
      unit :: Length Double
      unit = sq / (fromIntegral p *~ one)
      u = round ((d / unit) /~ one)
      (sqs, d1) = u `divMod` p
      s = show d1
      pad = if n == 0 then "" else replicate (n P.- length s) '0' ++ s
      