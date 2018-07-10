{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.Altitude (
  HasAltitude (..),
  groundPosition
) where

import Control.Lens(Lens', set)
import Numeric.Units.Dimensional.Prelude(Length, _0)

-- | All geographical coordinate systems need the concept of#
-- altitude above a reference point, usually associated with
-- local sea level.
-- 
-- Minimum definition: altitude, setAltitude.
class HasAltitude a where
  altitude ::
    Lens' a (Length Double)

instance HasAltitude (Length Double) where
  altitude =
    id

groundPosition ::
  HasAltitude a =>
  a
  -> a
groundPosition =
  set altitude _0
