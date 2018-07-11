{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.Longitude (
  HasLongitude (..),
) where

import Control.Lens(Lens')
import Data.Function(id)
import Numeric.Units.Dimensional.Prelude(Angle)
import Prelude(Double)

class HasLongitude a where
  longitudeL ::
    Lens' a (Angle Double)

instance HasLongitude (Angle Double) where
  longitudeL =
    id
