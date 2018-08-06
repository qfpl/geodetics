module Geodetics.Geodetic (
   -- ** Geodetic Coordinates
   Geodetic (..),
   HasGeodetic(..),
   readGroundPosition,
   toLocal,
   toWGS84,
   antipode,
   geometricalDistance,
   geometricalDistanceSq,
   groundDistance,
   properAngle,
   showAngle,
   -- ** Earth Centred Earth Fixed Coordinates
   ECEF,
   geoToEarth,
   earthToGeo,
   -- ** Re-exported for convenience
   _WGS84
) where

import Control.Lens(Lens', (^.), _Wrapped')
import Data.Char (chr)
import Data.Function((.))
import Data.Maybe
import Data.Monoid
import Geodetics.Ellipsoids
import Geodetics.LatLongParser
import Geodetics.Types.Altitude
import Geodetics.Types.Latitude
import Geodetics.Types.Longitude
import Geodetics.Types.Ellipsoid
import Geodetics.Types.TRF
import Linear.V3(V3(V3))
import Numeric.Units.Dimensional.Prelude hiding ((.))
import Text.ParserCombinators.ReadP
import qualified Prelude as P

-- | Defines a three-D position on or around the Earth using latitude,
-- longitude and altitude with respect to a specified ellipsoid, with
-- positive directions being North and East.  The default "show"
-- instance gives position in degrees, minutes and seconds to 5 decimal 
-- places, which is a
-- resolution of about 1m on the Earth's surface. Internally latitude
-- and longitude are stored as double precision radians. Convert to
-- degrees using e.g.  @latitude g /~ degree@.
-- 
-- The functions here deal with altitude by assuming that the local
-- height datum is always co-incident with the ellipsoid in use,
-- even though the \"mean sea level\" (the usual height datum) can be tens
-- of meters above or below the ellipsoid, and two ellipsoids can
-- differ by similar amounts. This is because the altitude is
-- usually known with reference to a local datum regardless of the
-- ellipsoid in use, so it is simpler to preserve the altitude across 
-- all operations. However if
-- you are working with ECEF coordinates from some other source then
-- this may give you the wrong results, depending on the altitude
-- correction your source has used.
-- 
-- There is no "Eq" instance because comparing two arbitrary
-- co-ordinates on the Earth is a non-trivial exercise. Clearly if all
-- the parameters are equal on the same ellipsoid then they are indeed
-- in the same place. However if different ellipsoids are used then
-- two co-ordinates with different numbers can still refer to the same
-- physical location.  If you want to find out if two co-ordinates are
-- the same to within a given tolerance then use "geometricDistance"
-- (or its squared variant to avoid an extra @sqrt@ operation).
data Geodetic =
  Geodetic
     Latitude
     Longitude
     Altitude
     TRF

class HasGeodetic a where
   geodetic ::
     Lens' a Geodetic

instance HasGeodetic Geodetic where
   geodetic =
      id

instance ManyTRF Geodetic where
  _ManyTRF =
    trf

instance FoldTRF Geodetic where
  _FoldTRF =
    trf

instance GetTRF Geodetic where
  _GetTRF =
    trf

instance SetTRF Geodetic where
  _SetTRF =
    trf

instance HasEllipsoid Geodetic where
  ellipsoid =
    trf . ellipsoid

instance GetEllipsoid Geodetic where
  _GetEllipsoid =
    ellipsoid

instance ManyEllipsoid Geodetic where
  _ManyEllipsoid =
    ellipsoid

instance FoldEllipsoid Geodetic where
  _FoldEllipsoid =
    ellipsoid

instance SetEllipsoid Geodetic where
  _SetEllipsoid =
    ellipsoid

instance HasTRF Geodetic where
   trf k (Geodetic lat' lon' alt' e) =
      fmap (\x -> Geodetic lat' lon' alt' x) (k e)

instance HasLatitude Geodetic where
   latitude k (Geodetic lat' lon' alt' e) =
      fmap (\x -> Geodetic x lon' alt' e) (k lat')

instance HasLongitude Geodetic where
   longitude k (Geodetic lat' lon' alt' e) =
      fmap (\x -> Geodetic lat' x alt' e) (k lon')

instance HasAltitude Geodetic where
  altitude k (Geodetic lat' lon' alt' e) =
      fmap (\x -> Geodetic lat' lon' x e) (k alt')

instance Show Geodetic where
   show g = concat [
      showAngle (abs $ (g ^. latitude . _Wrapped')),  " ", letter "SN" (g ^. latitude . _Wrapped'),  ", ",
      showAngle (abs $ (g ^. longitude . _Wrapped')), " ", letter "WE" (g ^. longitude . _Wrapped'), ", ", 
      show (g ^. altitude), " ", show (g ^. trf)]
      where letter s n = [s !! (if n < _0 then 0 else 1)]



-- | Read the latitude and longitude of a ground position and 
-- return a Geodetic position on the specified ellipsoid.
-- 
-- The latitude and longitude may be in any of the following formats.
-- The comma between latitude and longitude is optional in all cases. 
-- Latitude must always be first.
-- 
-- * Signed decimal degrees: 34.52327, -46.23234
-- 
-- * Decimal degrees NSEW: 34.52327N, 46.23234W
--
-- * Degrees and decimal minutes (units optional): 34° 31.43' N, 46° 13.92' 
-- 
-- * Degrees, minutes and seconds (units optional): 34° 31' 23.52\" N, 46° 13' 56.43\" W 
-- 
-- * DDDMMSS format with optional leading zeros: 343123.52N, 0461356.43W
readGroundPosition :: TRF -> String -> Maybe Geodetic
readGroundPosition e str = 
   case map fst $ filter (null . snd) $ readP_to_S latLong str of
      [] -> Nothing
      (lat,long) : _ -> Just $ groundPosition $ Geodetic (Latitude (lat *~ degree)) (Longitude (long *~ degree)) undefined e
      
      
-- | Show an angle as degrees, minutes and seconds to two decimal places.
showAngle :: Angle Double -> String
showAngle a
   | isNaN a1       = "NaN"  -- Not a Nangle
   | isInfinite a1  = sgn ++ "Infinity"
   | otherwise      = concat [sgn, show d, [chr 0xB0, ' '], 
                              show m, "' ", 
                              show s, ".", dstr, "\"" ]
   where
      a1 = a /~ one
      sgn = if a < _0 then "-" else ""
      centisecs :: Integer
      centisecs = P.abs $ P.round $ (a /~ degree) P.* 360000  -- hundredths of arcsec per degree.
      (d, m1) = centisecs `P.divMod` 360000
      (m, s1) = m1 `P.divMod` 6000   -- hundredths of arcsec per arcmin
      (s, ds) = s1 `P.divMod` 100
      dstr = reverse $ take 2 $ reverse (show ds) ++ "00" -- Decimal fraction with zero padding.
         
-- | The point on the Earth diametrically opposite the argument, with
-- the same altitude.
antipode :: Geodetic -> Geodetic
antipode g = Geodetic lat long (g ^. altitude) (g ^. trf)
   where
      lat = Latitude (negate $ (g ^. latitude . _Wrapped'))
      long' = (g ^. longitude . _Wrapped') - 180 *~ degree
      long | long' < _0  = Longitude (long' + 360 *~ degree)
           | otherwise  = Longitude long' 

   
   
-- | Convert a geodetic coordinate into earth centered, relative to the
-- ellipsoid in use.
geoToEarth :: Geodetic -> ECEF
geoToEarth geo = V3 (
      (n + h) * coslat * coslong)
      ((n + h) * coslat * sinlong)
      ((n * (_1 - eccentricity2 e) + h) * sinlat)
   where 
      geolat = geo ^. latitude . _Wrapped'
      geolon = geo ^. longitude . _Wrapped'
      n = normal e $ geolat
      e = geo ^. trf
      coslat = cos $ geolat
      coslong = cos $ geolon
      sinlat = sin $ geolat
      sinlong = sin $ geolon
      h = geo ^. altitude . _Wrapped'


-- | Convert an earth centred coordinate into a geodetic coordinate on 
-- the specified geoid.
--
-- Uses the closed form solution of H. Vermeille: Direct
-- transformation from geocentric coordinates to geodetic coordinates.
-- Journal of Geodesy Volume 76, Number 8 (2002), 451-454
earthToGeo :: TRF -> ECEF -> (Angle Double, Angle Double, Length Double)
earthToGeo e (V3 x y z) = (phi, atan2 y x, sqrt (l ^ pos2 + p2) - norm)
   where
      -- Naming: numeric suffix inicates power. Hence x2 = x * x, x3 = x2 * x, etc.
      p2 = x ^ pos2 + y ^ pos2
      a = (^. majorRadius) e
      a2 = a ^ pos2
      e2 = eccentricity2 e
      e4 = e2 ^ pos2
      zeta = (_1-e2) * (z ^ pos2 / a2)
      rho = (p2 / a2 + zeta - e4) / _6
      rho2 = rho ^ pos2
      rho3 = rho * rho2
      s = e4 * zeta * p2 / (_4 * a2)
      t = cbrt (s + rho3 + sqrt (s * (s + _2 * rho3)))
      u = rho + t + rho2 / t
      v = sqrt (u ^ pos2 + e4 * zeta)
      w = e2 * (u + v - zeta) / (_2 * v)
      kappa = _1 + e2 * (sqrt (u + v + w ^ pos2) + w) / (u + v)
      phi = atan (kappa * z / sqrt p2)
      norm = normal e phi
      l = z + e2 * norm * sin phi


-- | Convert a position from any geodetic to another one, assuming local altitude stays constant.
toLocal :: TRF -> Geodetic -> Geodetic
toLocal e2 g = Geodetic (Latitude lat) (Longitude lon) alt e2
   where
      alt = g ^. altitude
      (lat, lon, _) = earthToGeo e2 $ applyHelmert h $ geoToEarth g
      h = (^. helmert) (g ^. trf) `mappend` inverseHelmert ((^. helmert) e2)

-- | Convert a position from any geodetic to WGS84, assuming local
-- altitude stays constant.
toWGS84 :: Geodetic -> Geodetic
toWGS84 g = Geodetic (Latitude lat) (Longitude lon) alt _WGS84
   where
      alt = g ^. altitude
      (lat, lon, _) = earthToGeo _WGS84 $ applyHelmert h $ geoToEarth g
      h = (^. helmert) (g ^. trf)


-- | The absolute distance in a straight line between two geodetic 
-- points. They must be on the same ellipsoid.
-- Note that this is not the geodetic distance taken by following 
-- the curvature of the earth.
geometricalDistance :: Geodetic -> Geodetic -> Length Double
geometricalDistance g1 g2 = sqrt $ geometricalDistanceSq g1 g2

-- | The square of the absolute distance. Comes out as "Area" type of course.
geometricalDistanceSq :: Geodetic -> Geodetic -> Area Double
geometricalDistanceSq g1 g2 = (x1-x2) ^ pos2 + (y1-y2) ^ pos2 + (z1-z2) ^ pos2
   where
      (V3 x1 y1 z1) = geoToEarth g1
      (V3 x2 y2 z2) = geoToEarth g2


-- | The shortest ellipsoidal distance between two points on the
-- ground with reference to the same ellipsoid. Altitude is ignored.
--
-- The results are the distance between the points, the bearing of
-- the second point from the first, and (180 degrees - the bearing
-- of the first point from the second).
--
-- The algorithm can fail to converge where the arguments are near to
-- antipodal. In this case it returns @Nothing@.
--
-- Uses Vincenty's formula. \"Direct and inverse solutions of
-- geodesics on the ellipsoid with application of nested
-- equations\". T. Vincenty. Survey Review XXII 176, April
-- 1975. <http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf>
groundDistance :: Geodetic -> Geodetic ->
                  Maybe (Length Double, Dimensionless Double, Dimensionless Double)
groundDistance p1 p2 = do
     (_, (lambda, (cos2Alpha, delta, sinDelta, cosDelta, cos2DeltaM))) <-
       listToMaybe $ dropWhile converging $ take 100 $ zip lambdas $ tail lambdas
     let
       uSq = cos2Alpha * (a^pos2 - b^pos2) / b^pos2
       bigA = _1 + uSq/(16384*~one) * ((4096*~one) + uSq *
                                      (((-768)*~one) + uSq * ((320*~one)
                                                            - (175*~one)*uSq)))
       bigB = uSq/(1024*~one) * ((256*~one) +
                                 uSq * (((-128)*~one) +
                                        uSq * ((74*~one) - (47*~one)*uSq)))
       deltaDelta =
         bigB * sinDelta * (cos2DeltaM +
                             bigB/_4 * (cosDelta * (_2 * cos2DeltaM^pos2 - _1)
                                        - bigB/_6 * cos2DeltaM * (_4 * sinDelta^pos2 - _3)
                                          * (_4 * cos2DeltaM - _3)))
       s = b * bigA * (delta - deltaDelta)
       alpha1 = atan2(cosU2 * sin lambda) (cosU1 * sinU2 - sinU1 * cosU2 * cos lambda)
       alpha2 = atan2(cosU1 * sin lambda) (cosU1 * sinU2 * cos lambda - sinU1 * cosU2)
     return (s, alpha1, alpha2)
  where
    p1ellipsoid = p1 ^. trf
    f = flattening p1ellipsoid
    a = (^. majorRadius) p1ellipsoid
    b = minorRadius p1ellipsoid
    l = abs $ (p1 ^. longitude . _Wrapped') - (p2 ^. longitude . _Wrapped')
    u1 = atan ((_1-f) * tan (p1 ^. latitude . _Wrapped'))
    u2 = atan ((_1-f) * tan (p2 ^. latitude . _Wrapped'))
    sinU1 = sin u1
    cosU1 = cos u1
    sinU2 = sin u2
    cosU2 = cos u2
    
    nextLambda lambda = (lambda1, (cos2Alpha, delta, sinDelta, cosDelta, cos2DeltaM))
      where
        sinLambda = sin lambda
        cosLambda = cos lambda
        sinDelta = sqrt((cosU2 * sinLambda) ^ pos2 +
                        (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) ^ pos2)
        cosDelta = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda
        delta = atan2 sinDelta cosDelta
        sinAlpha = cosU1 * cosU2 * sinLambda / sinDelta
        cos2Alpha = _1 - sinAlpha ^ pos2
        cos2DeltaM = if cos2Alpha == _0
                     then _0
                     else cosDelta - _2 * sinU1 * sinU2 / cos2Alpha
        c = f/(16 *~ one) * cos2Alpha * (_4 + f * (_4 - _3 * cos2Alpha))
        lambda1 = l + (_1-c) * f * sinAlpha
                  * (delta + c * sinDelta
                     * (cos2DeltaM + c * cosDelta *(_2 * cos2DeltaM ^ pos2 - _1)))
    lambdas = iterate (nextLambda . fst) (l, undefined)
    converging ((l1,_),(l2,_)) = abs (l1 - l2) > (1e-14 *~ one)


-- | Add or subtract multiples of 2*pi so that for all @t@, @-pi < properAngle t < pi@.
properAngle :: Angle Double -> Angle Double
properAngle t 
   | r1 <= negate pi    = r1 + pi2
   | r1 > pi            = r1 - pi2
   | otherwise          = r1 
   where
      pf :: Double -> (Int, Double)
      pf = properFraction  -- Shut up GHC warning about defaulting to Integer.
      (_,r) = pf (t/pi2 /~ one)
      r1 = (r *~ one) * pi2
      pi2 = pi * _2

