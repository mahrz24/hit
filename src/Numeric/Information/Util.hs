{- |
Module      :  Numeric.Information.Util
Description :  Utility functions
Copyright   :  (c) Malte Harder
License     :  MIT

Maintainer  :  malte.harder@gmail.com
Stability   :  experimental
Portability :  portable

This module contains various helper functions for a more efficient
handling of probability distributions
-}

module Numeric.Information.Util (
    fromListWith'
  , fromListWithKey'
  ) where

import qualified Data.Map as Map
import qualified Data.List as List

fromListWith' :: Ord k => (a -> a -> a) -> [(k,a)] -> Map.Map k a
fromListWith' f
  = fromListWithKey' (\_ x y -> f x y)

fromListWithKey' :: Ord k => (k -> a -> a -> a) -> [(k,a)] -> Map.Map k a 
fromListWithKey' f xs
  = List.foldl' ins Map.empty xs
  where
    ins t (k,x) = Map.insertWithKey' f k x t

