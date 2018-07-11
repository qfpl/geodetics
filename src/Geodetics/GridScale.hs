{-# LANGUAGE TypeSynonymInstances, FlexibleInstances #-}

module Geodetics.GridScale (
  HasGridScale (..),
) where

import Control.Lens(Lens')
import Numeric.Units.Dimensional(Dimensionless)

class HasGridScale a where
  gridScale ::
    Lens' a (Dimensionless Double)

instance HasGridScale (Dimensionless Double) where
  gridScale =
    id
