{- |
Module      :  Numeric.Information.IT
Description :  Information quantities
Copyright   :  (c) Malte Harder
License     :  MIT

Maintainer  :  malte.harder@gmail.com
Stability   :  experimental
Portability :  portable

This module contains various function to calculate information
theoretic quantities
-}

module Numeric.Information.IT (
    -- * Expected values
    expected
    -- * Entropies
  , entropy
  , condEntropy
    -- * Divergences
  , divKL
  , divJS
    -- * Mutual Information
  , mutualInfo
  , mutualInfoJ
  , condMutualInfo
    -- * Information flow
    -- * Alternative Implementations
    -- ** Expected values
  , expected'
    -- ** Entropies
  , entropy'
  , condEntropy'
    -- * Mutual Information
  , mutualInfo'
  --  -- * Information flow
  --, informationFlow
    -- * Helper functions
  , log2
  , eventDivergence
  ) where

import Numeric.Information.Distribution
import Numeric.Probability.Distribution hiding (expected)

import qualified Data.List as List
import Control.Monad (liftM, liftM2)
import Control.Arrow (first, second)

log2 :: (Floating prob)
        => prob
        -> prob
log2 = logBase 2

minVal :: (Fractional prob) => prob
minVal = 1e-15

-- | Summing a function over the events of a distribution
-- weighted by the events probabilities
expected :: (Fractional prob, Ord prob, Ord a)
            => Dist prob a -- ^ p(x)
            -> (a -> prob) -- ^ f : X -> R
            -> prob -- ^ Sum_x p(x)*f(x)
expected p q = sum' . List.map (uncurry (*)) $
               List.filter ((> minVal) . snd) $
               (List.map . first $ q) $
               decons (normS p)

-- | Entropy of a distribution H(A)
entropy :: (Floating prob, Ord prob, Ord a)
           => Dist prob a
           -> prob
entropy p = - expected p (\a -> log2 $! a `seq` p ?= a)

-- | Conditional entropy H(A|B)
condEntropy :: (Floating prob, Ord prob, Ord a, Ord b)
               => (a -> Dist prob b)
               -> Dist prob a
               -> prob
condEntropy p_ygx p_x =
  let p_yx = p_ygx -|- p_x
  in - expected p_yx (\(y,x) -> log2 $ p_ygx x ?= y)

-- | Kullback Leibler Divergence of two distributions
divKL :: (Floating prob, Ord prob, Ord a)
         => Dist prob a
         -> Dist prob a
         -> prob
divKL p q = expected p (\a -> log2 $ (p ?= a)/(q ?= a))

-- | Jensen Shannon Divergence
divJS :: (Floating prob, Ord prob, Ord a)
         => prob -- ^ pi
         -> Dist prob a
         -> Dist prob a
         -> prob
divJS pi p q =
  let m = combine (convex pi) p q
  in entropy m - convex pi (entropy p) (entropy q)

-- | Mutual information of two random variables
mutualInfo :: (Floating prob, Ord prob, Ord a, Ord b)
              => (a -> Dist prob b) -- ^
              -> Dist prob a
              -> prob
mutualInfo p_ygx p_x =
  let p_yx = p_ygx -|- p_x
      p_y = margin fst p_yx
  in expected p_yx (\(y,x) ->
                     log2 $ eventDivergence p_yx p_y p_x (y,x))

-- | Mutual information of two random variables given by a joint distribution
mutualInfoJ :: (Floating prob, Ord prob, Ord a, Ord b)
               => Dist prob (a,b) -- ^
               -> prob
mutualInfoJ p_xy =
  let p_x = margin fst p_xy
      p_y = margin snd p_xy
  in expected p_xy (\(x,y) -> log2 $ eventDivergence p_xy p_x p_y (x,y))

-- | Conditional mutual information I(A;B|C)
condMutualInfo :: (Floating prob, Ord prob, Ord a, Ord b, Ord c)
                  => (c -> Dist prob a) -- ^
                  -> ((b,c) -> Dist prob a)
                  -> Dist prob c
                  -> Dist prob (b,c)
                  -> prob
condMutualInfo p_xgz p_xgyz p_z p_yz
  = condEntropy p_xgz p_z - condEntropy p_xgyz p_yz

{-
-- | Information flow between two random variables
informationFlow :: (Floating prob, Ord prob, Ord a, Ord b)
                   => (a -> Dist prob b) -- ^ Interventional distribution
                   -> Dist prob a
                   -> prob
informationFlow p_ygx p_x =
  let p_yx = p_ygx -|- p_x
  in expected p_yx (\(y,x) ->
                     log2 $ (p_ygx x ?= y)/
                     (expected p_x (\x' -> (p_ygx x' ?= y) )))

-- | Imposed information flow
imposedInformationFlow :: (Floating prob, Ord prob, Ord a, Ord b, Ord c)
                          => ((b,c) -> Dist prob a) -- ^ Interventional distribution
                          -> (c -> Dist prob b) -- ^ Interventional distribution
                          -> Dist prob c
                          -> prob
imposedInformationFlow p_xgyz p_ygz p_z =
  let p_yz = p_ygz -|- p_z
  in expected p_yz (\(y,z) -> expected (p_xgyz (y,z))
                              (\x ->
                                log2 $ (p_xgyz (y,z) ?= x)/
                                (expected (p_ygz z)
                                 (\y' -> (p_xgyz (y',z) ?= x)))))
-}

-- Alternative definitions of entropy and mutual information
-- maybe more efficient to use

-- | Summing a function over the events of a distribution,
-- weighted by the events probabilities with custom combinator
expected' :: (Fractional prob, Ord prob, Ord a)
             => ((prob,prob) -> prob) -- ^
             -> Dist prob a
             -> (a -> prob)
             -> prob
expected' f p q =
  sum' . List.map f $ List.filter ((> minVal) . snd) $
  (List.map . first $ q) $ decons (normS p)

entropy' :: (Floating prob, Ord prob, Ord a)
            => Dist prob a -- ^
            -> prob
entropy' p = - expected' (\(c,d) -> d * log2 d) p (const 1.0)

condEntropy' :: (Floating prob, Ord prob, Ord a, Ord b)
                => (a -> Dist prob b) -- ^
                -> Dist prob a
                -> prob
condEntropy' p_ygx p_x =
  let p_yx = p_ygx -|- p_x
  in - expected' (\(c,d) -> d * log2 (d*c)) p_yx (\(y,x) -> 1.0/(p_x ?= x) )

mutualInfo' :: (Floating prob, Ord prob, Ord a, Ord b)
               => (a -> Dist prob b)
               -> Dist prob a
               -> prob
mutualInfo' p_ygx p_x =
  let p_yx = p_ygx -|- p_x
      p_y = margin fst p_yx
  in expected' (\(c,d) -> d * log2 (d*c)) p_yx (\(y,x) -> 1.0/((p_x ?= x)*(p_y ?= y)))



-- Helper function
eventDivergence :: (Fractional prob, Ord a, Ord b)
                   => Dist prob (a,b)
                   -> Dist prob a
                   -> Dist prob b
                   -> (a,b)
                   -> prob
eventDivergence p_xy p_x p_y (x,y) =
  (p_xy ?= (x,y)) / ((p_x ?= x)*(p_y ?= y))
