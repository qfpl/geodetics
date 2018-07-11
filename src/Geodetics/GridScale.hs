{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.GridScale (
  HasGridScale (..),
) where

import Control.Lens(Lens')
import Data.Function(id)
import Numeric.Units.Dimensional(Dimensionless)
import Prelude(Double)

class HasGridScale a where
  gridScale ::
    Lens' a (Dimensionless Double)

instance HasGridScale (Dimensionless Double) where
  gridScale =
    id
