{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.Latitude (
  HasLatitude (..)
) where

import Control.Lens(Lens')
import Numeric.Units.Dimensional.Prelude(Angle)

class HasLatitude a where
  latitudeL ::
    Lens' a (Angle Double)

instance HasLatitude (Angle Double) where
  latitudeL =
    id
