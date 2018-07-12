{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE FlexibleContexts, TypeOperators, TypeFamilies #-}

{- | An TRF is a reasonable best fit for the surface of the 
Earth over some defined area. WGS84 is the standard used for the whole
of the Earth. Other Ellipsoids are considered a best fit for some
specific area.
-}

module Geodetics.Ellipsoids (
   -- ** Helmert transform between geodetic reference systems
   Helmert (..),
   HasHelmert(..),
   inverseHelmert,
   ECEF,
   applyHelmert,
   -- ** TRF models of the Geoid
   TRF (..),
   HasTRF(..),
   helmertFromWSG84,
   helmertToWSG84,
   _WGS84,
   flattening,
   minorRadius,
   eccentricity2,
   eccentricity'2,
   -- ** Auxiliary latitudes and related Values
   normal,
   latitudeRadius,
   meridianRadius,
   primeVerticalRadius,
   isometricLatitude,
   -- ** Tiny linear algebra library for 3D vectors
   V3,
   Matrix3,
   add3,
   mult3,
   scale3,
   negate3,
   transform3,
   invert3,
   trans3,
   dot3,
   cross3
) where

import Control.Lens(Lens', (^.))
import Control.Monad.Zip(MonadZip(mzipWith))
import Data.Monoid (Monoid)
import Data.Semigroup (Semigroup, (<>))
import Linear.V3(V3(V3))
import Numeric.Units.Dimensional.Prelude

-- | 3x3 transform matrix for V3.
type Matrix3 a = V3 (V3 a)

-- | Multiply a vector by a scalar.
scale3 :: (Num a, Functor f) =>
   f (Quantity d a) -> Quantity d' a -> f (Quantity (d * d') a)
scale3 v3 s = fmap (*s) v3

-- | Negation of a vector.
-- negate3 :: (Functor f, Num a) => f (Quantity d a) -> f (Quantity d a)
negate3 :: (Functor f, Num a) => f (Quantity d a) -> f (Quantity d a)
negate3 = fmap negate

-- | Add two vectors
add3 :: (Num a, MonadZip m) => m (Quantity d a) -> m (Quantity d a) -> m (Quantity d a)
add3 = mzipWith (+)

-- | Add two vectors
mult3 :: (Num a, MonadZip m) => m (Quantity d a) -> m (Quantity d a) -> m (Quantity (d * d) a)
mult3 = mzipWith (*)

-- | Multiply a matrix by a vector in the Dimensional type system.
transform3 :: (Num a) =>
   Matrix3 (Quantity d a) -> V3 (Quantity d' a) -> V3 (Quantity (d*d') a)
transform3 v3 v = fmap (\u -> dot3 u v) v3
   
-- | Inverse of a 3x3 matrix.
invert3 :: (Fractional a) =>
   Matrix3 (Quantity d a) -> Matrix3 (Quantity ((d*d)/(d*d*d)) a)
invert3 (V3 (V3 x1 y1 z1)
            (V3 x2 y2 z2)
            (V3 x3 y3 z3)) =
      V3
         (V3 (det2 y2 z2 y3 z3 / det) (det2 z1 y1 z3 y3 / det) (det2 y1 z1 y2 z2 / det))
         (V3 (det2 z2 x2 z3 x3 / det) (det2 x1 z1 x3 z3 / det) (det2 z1 x1 z2 x2 / det))
         (V3 (det2 x2 y2 x3 y3 / det) (det2 y1 x1 y3 x3 / det) (det2 x1 y1 x2 y2 / det))
   where
      det = (x1*y2*z3 + y1*z2*x3 + z1*x2*y3) - (z1*y2*x3 + y1*x2*z3 + x1*z2*y3)
      det2 a b c d = a*d - b*c

-- | Transpose of a 3x3 matrix.
trans3 :: Matrix3 a -> Matrix3 a
trans3 (V3 (V3 x1 y1 z1) (V3 x2 y2 z2) (V3 x3 y3 z3)) = V3 (V3 x1 x2 x3) (V3 y1 y2 y3) (V3 z1 z2 z3)


-- | Dot product of two vectors
dot3 :: (Num a) =>
   V3 (Quantity d1 a) -> V3 (Quantity d2 a) -> Quantity (d1 * d2) a
dot3 v3x v3y = sum (mzipWith (*) v3x v3y)

-- | Cross product of two vectors
cross3 :: (Num a) =>
   V3 (Quantity d1 a) -> V3 (Quantity d2 a) -> V3 (Quantity (d1 * d2) a)
cross3 (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (y1*z2 - z1*y2) (z1*x2 - x1*z2) (x1*y2 - y1*x2)


-- | The 7 parameter Helmert transformation. The monoid instance allows composition.
data Helmert = Helmert
   (Length Double)
   (Length Double)
   (Length Double)
   (Dimensionless Double)  -- Parts per million
   (Dimensionless Double)
   (Dimensionless Double)
   (Dimensionless Double)
   deriving (Eq, Show)

ccX ::
  HasHelmert a =>
  a
  -> Length Double
ccX =
   (^. cX)

ccY ::
  HasHelmert a =>
  a
  -> Length Double
ccY =
   (^. cY)

ccZ ::
  HasHelmert a =>
  a
  -> Length Double
ccZ =
   (^. cZ)

helmertScaleX ::
  HasHelmert a =>
  a
  -> Dimensionless Double
helmertScaleX =
   (^. helmertScale)

rrX ::
  HasHelmert a =>
  a
  -> Dimensionless Double
rrX =
   (^. rX)

rrY ::
  HasHelmert a =>
  a
  -> Dimensionless Double
rrY =
   (^. rY)

rrZ ::
  HasHelmert a =>
  a
  -> Dimensionless Double
rrZ =
   (^. rZ)

class HasHelmert a where
   helmert ::
     Lens' a Helmert
   cX ::
     Lens' a (Length Double)
   {-# INLINE cX #-}
   cY ::
     Lens' a (Length Double)
   {-# INLINE cY #-}
   cZ ::
     Lens' a (Length Double)
   {-# INLINE cZ #-}
   helmertScale ::
     Lens' a (Dimensionless Double)
   {-# INLINE helmertScale #-}
   rX ::
     Lens' a (Dimensionless Double)
   {-# INLINE rX #-}
   rY ::
     Lens' a (Dimensionless Double)
   {-# INLINE rY #-}
   rZ ::
     Lens' a (Dimensionless Double)
   {-# INLINE rZ #-}
   cX =
      helmert . cX
   cY =
      helmert . cY
   cZ =
      helmert . cZ
   helmertScale =
      helmert . helmertScale
   rX =
      helmert . rX
   rY =
      helmert . rY
   rZ =
      helmert . rZ

instance HasHelmert Helmert where
   {-# INLINE cX #-}
   {-# INLINE cY #-}
   {-# INLINE cZ #-}
   {-# INLINE helmertScale #-}
   {-# INLINE rX #-}
   {-# INLINE rY #-}
   {-# INLINE rZ #-}
   helmert =
      id
   cX k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert x cY' cZ' helmertScale' rX' rY' rZ') (k cX')
   cY k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert cX' x cZ' helmertScale' rX' rY' rZ') (k cY')
   cZ k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert cX' cY' x helmertScale' rX' rY' rZ') (k cZ')
   helmertScale k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert cX' cY' cZ' x rX' rY' rZ') (k helmertScale')
   rX k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert cX' cY' cZ' helmertScale' x rY' rZ') (k rX')
   rY k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert cX' cY' cZ' helmertScale' rX' x rZ') (k rY')
   rZ k (Helmert cX' cY' cZ' helmertScale' rX' rY' rZ') =
      fmap (\x -> Helmert cX' cY' cZ' helmertScale' rX' rY' x) (k rZ')

instance Semigroup Helmert where
    h1 <> h2 = Helmert (ccX h1 + ccX h2) (ccY h1 + ccY h2) (ccZ h1 + ccZ h2)
                       (helmertScaleX h1 + helmertScaleX h2)
                       (rrX h1 + rrX h2) (rrY h1 + rrY h2) (rrZ h1 + rrZ h2)

instance Monoid Helmert where
   mempty = Helmert (0 *~ meter) (0 *~ meter) (0 *~ meter) _0 _0 _0 _0
   mappend = (<>)

-- | The inverse of a Helmert transformation.
inverseHelmert :: HasHelmert e => e -> Helmert
inverseHelmert h = Helmert (negate $ ccX h) (negate $ ccY h) (negate $ ccZ h) 
                           (negate $ helmertScaleX h) 
                           (negate $ rrX h) (negate $ rrY h) (negate $ rrZ h)


-- | Earth-centred, Earth-fixed coordinates as a vector. The origin and axes are
-- not defined: use with caution.
type ECEF = V3 (Length Double)

-- | Apply a Helmert transformation to earth-centered coordinates.
applyHelmert:: HasHelmert e => e -> ECEF -> ECEF
applyHelmert h (V3 x y z) = V3 (
      ccX h + s * (                x - rrZ h * y + rrY h * z))
      (ccY h + s * (        rrZ h  * x +        y - rrX h * z))
      (ccZ h + s * (negate (rrY h) * x + rrX h * y +        z))
   where
      s = _1 + helmertScaleX h * (1e-6 *~ one)


-- | An TRF is defined by the major radius and the inverse flattening (which define its shape), 
-- and its Helmert transform relative to WGS84 (which defines its position and orientation).
--
-- The inclusion of the Helmert parameters relative to WGS84 actually make this a Terrestrial 
-- Reference Frame (TRF), but the term "Ellipsoid" will be used in this library for readability.
--
-- Minimum definition: @majorRadius@, @flatR@ & @helmert@.
-- 
-- Laws:
-- 
-- > helmertToWGS84 = applyHelmert . helmert
-- > helmertFromWGS84 e . helmertToWGS84 e = id
data TRF =
   TRF
      (Length Double)         -- majorRadius
      (Dimensionless Double)  -- flatR
      Helmert                 -- helmert
   deriving (Eq, Show)

class HasTRF a where
   trf ::   
      Lens' a TRF
   majorRadius ::
      Lens' a (Length Double)
   {-# INLINE majorRadius #-}
   flatR ::
      Lens' a (Dimensionless Double)
   {-# INLINE flatR #-}
   majorRadius =
      trf . majorRadius
   flatR =  
      trf . flatR

instance HasTRF TRF where
   {-# INLINE majorRadius #-}
   {-# INLINE flatR #-}
   trf =
      id
   majorRadius k (TRF r f h) =
      fmap (\r' -> TRF r' f h) (k r)
   flatR k (TRF r f h) =
      fmap (\f' -> TRF r f' h) (k f)

instance HasHelmert TRF where
   helmert k (TRF r f h) =
      fmap (\h' -> TRF r f h') (k h)

_WGS84 ::
  TRF
_WGS84 =
   TRF
      (6378137.0 *~ meter)
      (298.257223563 *~ one)
      mempty

helmertToWSG84 :: HasHelmert e => e -> ECEF -> ECEF
   -- ^ The Helmert transform that will convert a position wrt 
   -- this ellipsoid into a position wrt WGS84.
helmertToWSG84 e = applyHelmert (e ^. helmert)
helmertFromWSG84 :: HasHelmert e => e -> ECEF -> ECEF
   -- ^ And its inverse.
helmertFromWSG84 e = applyHelmert (inverseHelmert (e ^. helmert))

-- | The WGS84 geoid, major radius 6378137.0 meters, flattening = 1 / 298.257223563
-- as defined in \"Technical Manual DMA TM 8358.1 - Datums, Ellipsoids, Grids, and 
-- Grid Reference Systems\" at the National Geospatial-Intelligence Agency (NGA).
-- 
-- The WGS84 has a special place in this library as the standard TRF against
-- which all others are defined.

   
-- | Flattening (f) of an ellipsoid.
flattening :: HasTRF e => e -> Dimensionless Double
flattening e = _1 / (e ^. flatR)

-- | The minor radius of an ellipsoid.
minorRadius :: HasTRF e => e -> Length Double
minorRadius e = (e ^. majorRadius) * (_1 - flattening e)


-- | The eccentricity squared of an ellipsoid.
eccentricity2 :: HasTRF e => e -> Dimensionless Double
eccentricity2 e = _2 * f - (f * f) where f = flattening e

-- | The second eccentricity squared of an ellipsoid.
eccentricity'2 :: HasTRF e => e -> Dimensionless Double
eccentricity'2 e = (f * (_2 - f)) / (_1 - f * f) where f = flattening e


-- | Distance from the surface at the specified latitude to the 
-- axis of the Earth straight down. Also known as the radius of 
-- curvature in the prime vertical, and often denoted @N@.
normal :: HasTRF e => e -> Angle Double -> Length Double
normal e lat = (e ^. majorRadius) / sqrt (_1 - eccentricity2 e * sin lat ^ pos2)


-- | Radius of the circle of latitude: the distance from a point 
-- at that latitude to the axis of the Earth.
latitudeRadius :: HasTRF e => e -> Angle Double -> Length Double
latitudeRadius e lat = normal e lat * cos lat


-- | Radius of curvature in the meridian at the specified latitude. 
-- Often denoted @M@.
meridianRadius :: HasTRF e => e -> Angle Double -> Length Double
meridianRadius e lat = 
   (e ^. majorRadius ) * (_1 - eccentricity2 e)
   / sqrt ((_1 - eccentricity2 e * sin lat ^ pos2) ^ pos3)
   

-- | Radius of curvature of the ellipsoid perpendicular to the meridian at the specified latitude.
primeVerticalRadius :: HasTRF e => e -> Angle Double -> Length Double
primeVerticalRadius e lat =
   (e ^. majorRadius ) / sqrt (_1 - eccentricity2 e * sin lat ^ pos2)


-- | The isometric latitude. The isometric latitude is conventionally denoted by ψ 
-- (not to be confused with the geocentric latitude): it is used in the development 
-- of the ellipsoidal versions of the normal Mercator projection and the Transverse 
-- Mercator projection. The name "isometric" arises from the fact that at any point 
-- on the ellipsoid equal increments of ψ and longitude λ give rise to equal distance 
-- displacements along the meridians and parallels respectively.
isometricLatitude :: HasTRF e => e -> Angle Double -> Angle Double
isometricLatitude ellipse lat = atanh sinLat - e * atanh (e * sinLat)
   where
      sinLat = sin lat
      e = sqrt $ eccentricity2 ellipse
