{-# LANGUAGE NoImplicitPrelude #-}

module Geodetics.TransverseMercator(
   GridTM (..),
   HasGridTM(..),
   mkGridTM
) where

<<<<<<< HEAD
import Control.Lens((^.), _Wrapped')
import Data.Function
import Data.Monoid
import Geodetics.Ellipsoids
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.Types.Altitude
import Geodetics.Types.Latitude
import Geodetics.Types.Longitude
import Geodetics.Types.Ellipsoid
import Geodetics.Types.TRF
import Numeric.Units.Dimensional.Prelude hiding ((.))
=======
import Control.Lens(Lens', (^.))
import Data.Function((.), id)
import Geodetics.Altitude(HasAltitude(altitude))
import Geodetics.Ellipsoids
import Geodetics.Geodetic(HasGeodetic(geodetic), Geodetic(Geodetic))
import Geodetics.Grid(HasGridOffset(gridOffsetL), GridClass(gridTRF), GridOffset(GridOffset), GridPoint(GridPoint), fromGrid, toGrid, applyOffset, gridBasis, offsetNegate)
import Geodetics.GridScale(HasGridScale(gridScale))
import Geodetics.Latitude(HasLatitude(latitude))
import Geodetics.Longitude(HasLongitude(longitude))
import Numeric.Units.Dimensional.Prelude hiding ((.), id)
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
import Prelude ()

-- | A Transverse Mercator projection gives an approximate mapping of the ellipsoid on to a 2-D grid. It models
-- a sheet curved around the ellipsoid so that it touches it at one north-south line (hence making it part of
-- a slightly elliptical cylinder).
data GridTM = GridTM
   Geodetic
      -- A point on the line where the projection touches the ellipsoid (altitude is ignored).
   GridOffset
      -- The grid position of the true origin. Used to avoid negative coordinates over 
      -- the area of interest. The altitude gives a vertical offset from the ellipsoid.
   (Dimensionless Double)
      -- A scaling factor that balances the distortion between the east & west edges and the middle 
      -- of the projection.
      
   -- Remaining elements are memoised parameters computed from the ellipsoid underlying the true origin.
   (Dimensionless Double)
   (Dimensionless Double)
   (Dimensionless Double)
   (Dimensionless Double)
   deriving (Show)

class HasGridTM a where
  gridTM ::
    Lens' a GridTM
  gridN1 ::
    Lens' a (Dimensionless Double)
  gridN2 ::
    Lens' a (Dimensionless Double)
  gridN3 ::
    Lens' a (Dimensionless Double)
  gridN4 ::
    Lens' a (Dimensionless Double)
  gridN1 =
    gridTM . gridN1
  {-# INLINE gridN1 #-}
  gridN2 =
    gridTM . gridN2
  {-# INLINE gridN2 #-}
  gridN3 =
    gridTM . gridN3
  {-# INLINE gridN3 #-}
  gridN4 =
    gridTM . gridN4
  {-# INLINE gridN4 #-}

instance HasGridTM GridTM where
  gridTM =
    id
  gridN1 k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM t f s x n2 n3 n4) (k n1)
  {-# INLINE gridN1 #-}
  gridN2 k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM t f s n1 x n3 n4) (k n2)
  {-# INLINE gridN2 #-}
  gridN3 k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM t f s n1 n2 x n4) (k n3)
  {-# INLINE gridN3 #-}
  gridN4 k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM t f s n1 n2 n3 x) (k n4)
  {-# INLINE gridN4 #-}

instance HasGridScale GridTM where
  gridScale k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM t f x n1 n2 n3 n4) (k s)

instance HasGeodetic GridTM where
  geodetic k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM x f s n1 n2 n3 n4) (k t)

instance HasGridOffset GridTM where
  gridOffsetL k (GridTM t f s n1 n2 n3 n4) =
    fmap (\x -> GridTM t x s n1 n2 n3 n4) (k f)

instance HasLatitude GridTM where
  latitude =
    geodetic . latitude

instance HasLongitude GridTM where
  longitude =
    geodetic . longitude

instance HasAltitude GridTM where
  altitude =
    geodetic . altitude

instance HasEllipsoid GridTM where
  ellipsoid =
    geodetic . ellipsoid

instance HasTRF GridTM where
  trf =
    geodetic . trf

-- | Create a Transverse Mercator grid.
mkGridTM :: 
   Geodetic                 -- ^ True origin.
   -> GridOffset            -- ^ Vector from true origin to false origin.
   -> Dimensionless Double  -- ^ Scale factor.
   -> GridTM
mkGridTM origin offset sf =
   GridTM  origin
           offset
           sf
           (_1 + n + (_5/_4) * n^pos2 + (_5/_4) * n^pos3)
           (_3 * n + _3 * n^pos2 + ((21*~one)/_8) * n^pos3)
           (((15*~one)/_8) * (n^pos2 + n^pos3))
           (((35*~one)/(24*~one)) * n^pos3)
    where 
       f = flattening (origin ^. trf)
       n = f / (_2-f)  -- Equivalent to (a-b)/(a+b) where b = (1-f)*a




-- | Equation C3 from reference [1].
m :: GridTM -> Dimensionless Double -> Length Double
m grid lat = bF0 * (grid ^. gridN1 * dLat 
                    - grid ^. gridN2 * sin dLat * cos sLat
                    + grid ^. gridN3 * sin (_2 * dLat) * cos (_2 * sLat) 
                    - grid ^. gridN4 * sin (_3 * dLat) * cos (_3 * sLat))
   where
<<<<<<< HEAD
      dLat = lat - (trueOrigin grid ^. latitude . _Wrapped')
      sLat = lat + (trueOrigin grid ^. latitude . _Wrapped')
      bF0 = minorRadius (gridEllipsoid grid) * gridScale grid
=======
      dLat = lat - (grid ^. geodetic . latitude)
      sLat = lat + (grid ^. geodetic . latitude)
      bF0 = minorRadius (gridTRF grid) * (grid ^. gridScale)
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d


instance GridClass GridTM where
   fromGrid p = Geodetic
      (Latitude (lat' - east' ^ pos2 * tanLat / (_2 * rho * v)  -- Term VII
            + east' ^ pos4 * (tanLat / ((24 *~ one) * rho * v ^ pos3)) 
                           * (_5 + _3 * tanLat ^ pos2 + eta2 - _9 * tanLat ^ pos2 * eta2)  -- Term VIII
            - east' * east' ^ pos5 * (tanLat / ((720 *~ one) * rho * v ^ pos5))
<<<<<<< HEAD
                           * (61 *~ one + (90 *~ one) * tanLat ^ pos2 + (45 *~ one) * tanLat ^ pos4))) -- Term IX
      (Longitude ((trueOrigin grid ^. longitude . _Wrapped') 
=======
                           * (61 *~ one + (90 *~ one) * tanLat ^ pos2 + (45 *~ one) * tanLat ^ pos4)) -- Term IX
      ((grid ^. geodetic . longitude) 
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
            + east' / (cosLat * v)  -- Term X
            - (east' ^ pos3 / (_6 * cosLat * v ^ pos3)) * (v / rho + _2 * tanLat ^ pos2)  -- Term XI
            + (east' ^ pos5 / ((120 *~ one) * cosLat * v ^ pos5)) 
                 * (_5 + (28 *~ one) * tanLat ^ pos2  + (24 *~ one) * tanLat ^ pos4)  -- Term XII
            - (east' ^ pos5 * east' ^ pos2 / ((5040 *~ one) * cosLat * v * v * v ^ pos5))
<<<<<<< HEAD
                 * ((61 *~ one) + (662 *~ one) * tanLat ^ pos2 + (1320 *~ one) * tanLat ^ pos4 + (720 *~ one) * tanLat * tanLat ^ pos5))) -- Term XIIa
     (Altitude (0 *~ meter)) (gridEllipsoid grid)
=======
                 * ((61 *~ one) + (662 *~ one) * tanLat ^ pos2 + (1320 *~ one) * tanLat ^ pos4 + (720 *~ one) * tanLat * tanLat ^ pos5)) -- Term XIIa
     (0 *~ meter) (gridTRF grid)
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
            
            
      where
         GridPoint east' north' _ _ = (grid ^. gridOffsetL) `applyOffset` p
         lat' = fst $ head $ dropWhile ((> 0.01 *~ milli meter) . snd) 
<<<<<<< HEAD
               $ tail $ iterate next (trueOrigin grid ^. latitude . _Wrapped', 1 *~ meter) 
=======
               $ tail $ iterate next (grid ^. geodetic . latitude, 1 *~ meter) 
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
            where
               next (phi, _) = let delta = north' - m grid phi in (phi + delta / aF0, delta) 
               -- head and tail are safe because iterate returns an infinite list.
          
         sinLat = sin lat'
         cosLat = cos lat'
         tanLat = tan lat'
         sinLat2 = sinLat ^ pos2
         v = aF0 / sqrt (_1 - e2 * sinLat2)
         rho = aF0 * (_1 - e2) * (_1 - e2 * sinLat2) ** ((-1.5) *~ one)
         eta2 = v / rho - _1
               
               
<<<<<<< HEAD
         aF0 = (^. majorRadius) (gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         grid = gridBasis p
=======
         aF0 = (gridTRF grid ^. majorRadius ) * (grid ^. gridScale)
         e2 = eccentricity2 $ gridTRF grid
         grid = p ^. gridBasis
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
         
   toGrid grid geo = applyOffset (off  `mappend` (offsetNegate (grid ^. gridOffsetL))) $ 
                     GridPoint _0 _0 _0 grid
      where
         v = aF0 / sqrt (_1 - e2 * sinLat2)
         rho = aF0 * (_1 - e2) * (_1 - e2 * sinLat2) ** ((-1.5) *~ one)
         eta2 = v / rho - _1
         off = GridOffset
                  (dLong * term_IV
                   + dLong ^ pos3 * term_V
                   + dLong ^ pos5 * term_VI)
                  (m grid lat + dLong ^ pos2 * term_II
                     + dLong ^ pos4 * term_III 
                     + dLong * dLong ^ pos5 * term_IIIa)
                  (0 *~ meter)
         -- Terms defined in [1].
         term_II   = (v/_2) * sinLat * cosLat
         term_III  = (v/(24*~one)) * sinLat * cosLat ^ pos3 
                     * (_5 - tanLat ^ pos2 + _9 * eta2)
         term_IIIa = (v/(720*~one)) * sinLat * cosLat ^ pos5 
                     * ((61 *~ one) - (58 *~ one) * tanLat ^ pos2 + tanLat ^ pos4)
         term_IV   = v * cosLat
         term_V    = (v/_6) * cosLat ^ pos3 * (v/rho - tanLat ^ pos2)
         term_VI   = (v/(120*~one)) * cosLat ^ pos5 
                     * (_5 - (18*~one) * tanLat ^ pos2 
                              + tanLat ^ pos4 + (14*~one) * eta2
                              - (58*~one) * tanLat ^ pos2 * eta2)
         {- 
         -- Trace message for debugging. Uncomment this code for easy access to intermediate values.
         traceMsg = concat [
            "v    = ", show v, "\n",
            "rho  = ", show rho, "\n",
            "eta2 = ", show eta2, "\n",
            "M    = ", show $ m grid lat, "\n",
            "I    = ", show $ m grid lat + deltaNorth (falseOrigin grid), "\n",
            "II   = ", show term_II, "\n",
            "III  = ", show term_III, "\n",
            "IIIa = ", show term_IIIa, "\n",
            "IV   = ", show term_IV, "\n",
            "V    = ", show term_V, "\n",
            "VI   = ", show term_VI, "\n"]
         -}
         -- Common subexpressions
<<<<<<< HEAD
         lat = geo ^. latitude . _Wrapped'
         long = geo ^. longitude . _Wrapped'
         dLong = long - (trueOrigin grid ^. longitude . _Wrapped')
=======
         lat = geo ^. latitude
         long = geo ^. longitude
         dLong = long - (grid ^. geodetic . longitude)
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
         sinLat = sin lat
         cosLat = cos lat
         tanLat = tan lat
         sinLat2 = sinLat ^ pos2
<<<<<<< HEAD
         aF0 = ((^. majorRadius) $ gridEllipsoid grid) * gridScale grid
         e2 = eccentricity2 $ gridEllipsoid grid
         
   gridEllipsoid = (^. trf) . trueOrigin
=======
         aF0 = (gridTRF grid ^. majorRadius) * (grid ^. gridScale)
         e2 = eccentricity2 $ gridTRF grid
   
   gridTRF = (^. geodetic . trf)
>>>>>>> c5162f79af10ba19aeda66fac966dc31768d576d
