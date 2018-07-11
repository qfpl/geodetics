{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.Latitude (
  HasLatitude (..)
) where

import Control.Lens(Lens')
import Data.Function(id)
import Numeric.Units.Dimensional.Prelude(Angle)
import Prelude(Double)

class HasLatitude a where
  latitudeL ::
    Lens' a (Angle Double)

instance HasLatitude (Angle Double) where
  latitudeL =
    id
