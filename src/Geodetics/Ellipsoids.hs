{-# LANGUAGE FlexibleContexts, TypeOperators, TypeFamilies #-}

{- | An Ellipsoid is a reasonable best fit for the surface of the 
Earth over some defined area. WGS84 is the standard used for the whole
of the Earth. Other Ellipsoids are considered a best fit for some
specific area.
-}

module Geodetics.Ellipsoids (
   -- ** Helmert transform between geodetic reference systems
   Helmert (..),
   inverseHelmert,
   ECEF,
   applyHelmert,
   -- ** Ellipsoid models of the Geoid
   Ellipsoid (..),
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
   scale3,
   negate3,
   transform3,
   invert3,
   trans3,
   dot3,
   cross3
) where

import Data.Monoid (Monoid)
import Data.Semigroup (Semigroup, (<>))
import Linear.V3(V3(V3))
import Numeric.Units.Dimensional
import Numeric.Units.Dimensional.Prelude
import Prelude ()  -- Numeric instances.


-- | 3x3 transform matrix for V3.
type Matrix3 a = V3 (V3 a)

-- | Multiply a vector by a scalar.
scale3 :: (Num a) =>
   V3 (Quantity d a) -> Quantity d' a -> V3 (Quantity (d * d') a)
scale3 (V3 x y z) s = V3 (x*s) (y*s) (z*s)


-- | Negation of a vector.
negate3 :: (Num a) => V3 (Quantity d a) -> V3 (Quantity d a)
negate3 (V3 x y z) = V3 (negate x) (negate y) (negate z)

-- | Add two vectors
add3 :: (Num a) => V3 (Quantity d a) -> V3 (Quantity d a) -> V3 (Quantity d a)
add3 (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (x1+x2) (y1+y2) (z1+z2)


-- | Multiply a matrix by a vector in the Dimensional type system.
transform3 :: (Num a) =>
   Matrix3 (Quantity d a) -> V3 (Quantity d' a) -> V3 (Quantity (d*d') a)
transform3 (V3 tx ty tz) v = V3 (t tx v) (t ty v) (t tz v)
   where
      t (V3 x1 y1 z1) (V3 x2 y2 z2) = x1*x2 + y1*y2 + z1*z2


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
dot3 (V3 x1 y1 z1) (V3 x2 y2 z2) = x1*x2 + y1*y2 + z1*z2

-- | Cross product of two vectors
cross3 :: (Num a) =>
   V3 (Quantity d1 a) -> V3 (Quantity d2 a) -> V3 (Quantity (d1 * d2) a)
cross3 (V3 x1 y1 z1) (V3 x2 y2 z2) = V3 (y1*z2 - z1*y2) (z1*x2 - x1*z2) (x1*y2 - y1*x2)


-- | The 7 parameter Helmert transformation. The monoid instance allows composition.
data Helmert = Helmert {
   cX, cY, cZ :: Length Double,
   helmertScale :: Dimensionless Double,  -- ^ Parts per million
   rX, rY, rZ :: Dimensionless Double } deriving (Eq, Show)

instance Semigroup Helmert where
    h1 <> h2 = Helmert (cX h1 + cX h2) (cY h1 + cY h2) (cZ h1 + cZ h2)
                       (helmertScale h1 + helmertScale h2)
                       (rX h1 + rX h2) (rY h1 + rY h2) (rZ h1 + rZ h2)

instance Monoid Helmert where
   mempty = Helmert (0 *~ meter) (0 *~ meter) (0 *~ meter) _0 _0 _0 _0
   mappend = (<>)

-- | The inverse of a Helmert transformation.
inverseHelmert :: Helmert -> Helmert
inverseHelmert h = Helmert (negate $ cX h) (negate $ cY h) (negate $ cZ h) 
                           (negate $ helmertScale h) 
                           (negate $ rX h) (negate $ rY h) (negate $ rZ h)


-- | Earth-centred, Earth-fixed coordinates as a vector. The origin and axes are
-- not defined: use with caution.
type ECEF = V3 (Length Double)

-- | Apply a Helmert transformation to earth-centered coordinates.
applyHelmert:: Helmert -> ECEF -> ECEF
applyHelmert h (V3 x y z) = V3 (
      cX h + s * (                x - rZ h * y + rY h * z))
      (cY h + s * (        rZ h  * x +        y - rX h * z))
      (cZ h + s * (negate (rY h) * x + rX h * y +        z))
   where
      s = _1 + helmertScale h * (1e-6 *~ one)


-- | An Ellipsoid is defined by the major radius and the inverse flattening (which define its shape), 
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
data Ellipsoid =
   Ellipsoid {
      majorRadius ::
         Length Double
   ,  flatR ::
         Dimensionless Double
   ,  helmert ::
         Helmert
   }
   deriving (Eq, Show)

_WGS84 ::
  Ellipsoid
_WGS84 =
   Ellipsoid
      (6378137.0 *~ meter)
      (298.257223563 *~ one)
      mempty

helmertToWSG84 :: Ellipsoid -> ECEF -> ECEF
   -- ^ The Helmert transform that will convert a position wrt 
   -- this ellipsoid into a position wrt WGS84.
helmertToWSG84 e = applyHelmert (helmert e)
helmertFromWSG84 :: Ellipsoid -> ECEF -> ECEF
   -- ^ And its inverse.
helmertFromWSG84 e = applyHelmert (inverseHelmert $ helmert e)

-- | The WGS84 geoid, major radius 6378137.0 meters, flattening = 1 / 298.257223563
-- as defined in \"Technical Manual DMA TM 8358.1 - Datums, Ellipsoids, Grids, and 
-- Grid Reference Systems\" at the National Geospatial-Intelligence Agency (NGA).
-- 
-- The WGS84 has a special place in this library as the standard Ellipsoid against
-- which all others are defined.

   
-- | Flattening (f) of an ellipsoid.
flattening :: Ellipsoid -> Dimensionless Double
flattening e = _1 / flatR e

-- | The minor radius of an ellipsoid.
minorRadius :: Ellipsoid -> Length Double
minorRadius e = majorRadius e * (_1 - flattening e)


-- | The eccentricity squared of an ellipsoid.
eccentricity2 :: Ellipsoid -> Dimensionless Double
eccentricity2 e = _2 * f - (f * f) where f = flattening e

-- | The second eccentricity squared of an ellipsoid.
eccentricity'2 :: Ellipsoid -> Dimensionless Double
eccentricity'2 e = (f * (_2 - f)) / (_1 - f * f) where f = flattening e


-- | Distance from the surface at the specified latitude to the 
-- axis of the Earth straight down. Also known as the radius of 
-- curvature in the prime vertical, and often denoted @N@.
normal :: Ellipsoid -> Angle Double -> Length Double
normal e lat = majorRadius e / sqrt (_1 - eccentricity2 e * sin lat ^ pos2)


-- | Radius of the circle of latitude: the distance from a point 
-- at that latitude to the axis of the Earth.
latitudeRadius :: Ellipsoid -> Angle Double -> Length Double
latitudeRadius e lat = normal e lat * cos lat


-- | Radius of curvature in the meridian at the specified latitude. 
-- Often denoted @M@.
meridianRadius :: Ellipsoid -> Angle Double -> Length Double
meridianRadius e lat = 
   majorRadius e * (_1 - eccentricity2 e) 
   / sqrt ((_1 - eccentricity2 e * sin lat ^ pos2) ^ pos3)
   

-- | Radius of curvature of the ellipsoid perpendicular to the meridian at the specified latitude.
primeVerticalRadius :: Ellipsoid -> Angle Double -> Length Double
primeVerticalRadius e lat =
   majorRadius e / sqrt (_1 - eccentricity2 e * sin lat ^ pos2)


-- | The isometric latitude. The isometric latitude is conventionally denoted by ψ 
-- (not to be confused with the geocentric latitude): it is used in the development 
-- of the ellipsoidal versions of the normal Mercator projection and the Transverse 
-- Mercator projection. The name "isometric" arises from the fact that at any point 
-- on the ellipsoid equal increments of ψ and longitude λ give rise to equal distance 
-- displacements along the meridians and parallels respectively.
isometricLatitude :: Ellipsoid -> Angle Double -> Angle Double
isometricLatitude ellipse lat = atanh sinLat - e * atanh (e * sinLat)
   where
      sinLat = sin lat
      e = sqrt $ eccentricity2 ellipse
