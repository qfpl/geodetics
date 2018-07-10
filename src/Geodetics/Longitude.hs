{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.Longitude (
  HasLongitude (..),
) where

import Control.Lens(Lens')
import Numeric.Units.Dimensional.Prelude(Angle)

class HasLongitude a where
  longitudeL ::
    Lens' a (Angle Double)

instance HasLongitude (Angle Double) where
  longitudeL =
    id
