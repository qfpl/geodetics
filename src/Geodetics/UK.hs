{-# LANGUAGE MultiParamTypeClasses #-}

-- | Distinguished coordinate systems for the United Kingdom.
module Geodetics.UK (
   _OSGB36,
   UkNationalGrid (..),
   ukGrid,
   fromUkGridReference,
   toUkGridReference
) where

import Control.Applicative
import Control.Monad
import Data.Array
import Data.Char
import Data.Monoid
import Geodetics.Geodetic
import Geodetics.Grid
import Geodetics.TransverseMercator
import Geodetics.Types.Altitude
import Geodetics.Types.Latitude
import Geodetics.Types.Longitude
import Geodetics.Types.TRF
import Numeric.Units.Dimensional.Prelude
import qualified Prelude as P


-- | The UK National Grid is a Transverse Mercator projection with a true origin at
-- 49 degrees North, 2 degrees West on OSGB36, and a false origin 400km West and 100 km North of
-- the true origin. The scale factor is defined as @10**(0.9998268 - 1)@.
data UkNationalGrid = UkNationalGrid deriving (Eq, Show)

instance GridClass UkNationalGrid where
   toGrid _ = unsafeGridCoerce UkNationalGrid . toGrid ukGrid
   fromGrid = fromGrid . unsafeGridCoerce ukGrid
   gridEllipsoid _ = _OSGB36


ukTrueOrigin :: Geodetic
ukTrueOrigin =
  Geodetic
    (Latitude (49 *~ degree))
    (Longitude ((-2) *~ degree))
    (Altitude (0 *~ meter))
    (_OSGB36)

ukFalseOrigin :: GridOffset 
ukFalseOrigin = GridOffset ((-400) *~ kilo meter) (100 *~ kilo meter) (0 *~ meter)


-- | Numerical definition of the UK national grid.
ukGrid :: GridTM
ukGrid = mkGridTM ukTrueOrigin ukFalseOrigin 
   ((10 *~ one) ** (0.9998268 *~ one - _1))


-- | Size of a UK letter-pair grid square.
ukGridSquare :: Length Double
ukGridSquare = 100 *~ kilo meter


-- | Convert a grid reference to a position, if the reference is valid. 
-- This returns the position of the south-west corner of the nominated 
-- grid square and an offset to its centre. Altitude is set to zero.
fromUkGridReference :: String -> Maybe (GridPoint UkNationalGrid, GridOffset)
fromUkGridReference str = if length str < 2 then Nothing else do
      let 
         c1:c2:ds = str
         n = length ds
      guard $ even n
      let (dsE, dsN) = splitAt (n `div` 2) ds
      (east, sq) <- fromGridDigits ukGridSquare dsE
      (north, _) <- fromGridDigits ukGridSquare dsN
      base <- fromUkGridLetters c1 c2
      let half = sq / (2 *~ one)
      return (applyOffset (GridOffset east north (0 *~ meter)) base,
              GridOffset half half (0 *~ meter))

      


-- | The south west corner of the nominated grid square, if it is a legal square.
-- This function works for all pairs of letters except 'I' (as that is not used).
-- In practice only those pairs covering the UK are actually considered meaningful.
fromUkGridLetters :: Char -> Char -> Maybe (GridPoint UkNationalGrid)
fromUkGridLetters c1 c2 = applyOffset <$> (mappend <$> g1 <*> g2) <*> letterOrigin
   where
      letterOrigin = Just $ GridPoint ((-1000) *~ kilo meter) ((-500) *~ kilo meter) m0 UkNationalGrid
      gridIndex c = 
         if inRange ('A', 'H') c then Just $ ord c P.- ord 'A'  -- 'I' is not used.
         else if inRange ('J', 'Z') c then Just $ ord c P.- ord 'B'
         else Nothing
      gridSquare c = do -- Maybe monad
         g <- gridIndex c
         let (y,x) = g `divMod` 5 
         return (fromIntegral x *~ one, _4 - fromIntegral y *~ one)
      g1 = do
         (x,y) <- gridSquare c1
         return $ GridOffset (x * (500 *~ kilo meter)) (y * (500 *~ kilo meter)) m0
      g2 = do
         (x,y) <- gridSquare c2
         return $ GridOffset (x * (100 *~ kilo meter)) (y * (100 *~ kilo meter)) m0
      m0 = 0 *~ meter


-- | Find the nearest UK grid reference point to a specified position. The Int argument is the number of
-- digits precision, so 2 for a 4-figure reference and 3 for a 6-figure reference, although any value 
-- between 0 and 5 can be used (giving a 1 meter precision).
-- Altitude is ignored. If the result is outside the area defined by the two letter grid codes then
-- @Nothing@ is returned.
toUkGridReference :: Int -> GridPoint UkNationalGrid -> Maybe String
toUkGridReference n p
   | n < 0         = error "toUkGridReference: precision argument must not be negative."
   | otherwise     = do
      (gx, strEast) <- toGridDigits ukGridSquare n $ eastings p + 1000 *~ kilo meter
      (gy, strNorth) <- toGridDigits ukGridSquare n $ northings p + 500 *~ kilo meter
      let (gx1, gx2) = (fromIntegral gx) `divMod` 5
          (gy1, gy2) = (fromIntegral gy) `divMod` 5
      guard (gx1 < 5 && gy1 < 5)
      let c1 = gridSquare gx1 gy1
          c2 = gridSquare gx2 gy2
      return $ c1 : c2 : strEast ++ strNorth
   where
      gridSquare x y = letters ! (4 P.- y, x)
      letters :: Array (Int, Int) Char
      letters = listArray ((0,0),(4,4)) $ ['A'..'H'] ++ ['J'..'Z']
   