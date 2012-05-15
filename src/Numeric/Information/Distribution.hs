{-#LANGUAGE TypeFamilies #-}
{- |
Module      :  Numeric.Information.Distribution
Description :  Distribution calculus
Copyright   :  (c) Malte Harder
License     :  MIT

Maintainer  :  malte.harder@gmail.com
Stability   :  experimental
Portability :  portable

This module contains a few extensions to Numeric.Probability.Distribution
-}

module Numeric.Information.Distribution (
    -- * Data types
    Dist
    -- * Stricter versions of some Numeric.Probability.Distribution methods
  , uniform'
  , shape'
  , fromFreqs'
  , sumP'
  , sum'
    -- * Sorted normalization
  , normS
  , normS'
    -- * Distribution generators
  , delta
  , delta'
  , xor
    -- * Distribution accessor
  , (???)
  , (?=)
    -- * Distribution combinators
  , (<*>)
  , (-|-)
  , (<\/>)
  , (<=>)
  , (-<-)
    -- * Array index operators
  , (<:>)
  , (<**>)
  , (-|*-)
  , (<\/*>)
  , (<=*>)
  , (<=:>)
  , (-<*-)
    -- * Index independent combinators
  , (-\-)
  , (-.-)
  , (-:-)
  , combine
    -- * Marginalization
  , margin
    -- * Helper functions
  , convex
  , arrayIndex
    -- * Parsing
  , readDistribution
  , readCondDistribution
  ) where

import Numeric.Information.Util

import qualified Data.List as List
import qualified Data.Map as Map
import Numeric.Probability.Distribution as Prob
import Numeric.Probability.Shape as Shape
import Control.Monad (liftM,liftM2)
import Data.Maybe

-- Some extensions to the probability library

-- Distribution type
type Dist prob a = Prob.T prob a

uniform' :: Fractional prob
            => Spread prob a
uniform' = shape' Shape.uniform

shape' :: Fractional prob
          => (prob -> prob)
          -> Spread prob a
shape' _ [] = error "Probability.shape: empty list"
shape' f xs =
   let incr = 1 / fromIntegral (length xs - 1)
       ps = List.map f (iterate (+incr) 0)
   in  fromFreqs' (zip xs ps)

fromFreqs' :: (Fractional prob)
              => [(a,prob)]
              -> Dist prob a
fromFreqs' xs = Cons (List.map (\(x,p)->(x,p/q)) xs)
           where q = sumP' xs

sumP' :: Num prob
         => [(a,prob)]
         -> prob
sumP' = sum' . List.map snd

sum' :: Num a => [a] -> a
sum' = List.foldl' (+) 0

-- | Normalize and sort
normS :: (Num prob, Ord a)
         => Dist prob a
         -> Dist prob a
normS = lift normS'

normS' :: (Num prob, Ord a)
          => [(a,prob)]
          -> [(a,prob)]
normS' =
   Map.toAscList . fromListWith' (+)

delta :: (Fractional prob)
         => a -> Dist prob a
delta i = fromFreqs' [(i,1)]

delta' :: (Fractional prob)
         => [a] -> Dist prob a
delta' i = fromFreqs' [(head i,1)]

xor :: (Fractional prob)
       => (Bool, Bool)
       -> Dist prob Bool
xor (a,b) = fromFreqs' [((a/=b)==True,1)]



-- | Extract probability of event
(???) :: Num prob
         => Event a
         -> Dist prob a -- ^
         -> prob
(???) p = sumP' . List.filter (p . fst) . decons

-- | The probability of a certain event
(?=) :: (Num prob, Eq a)
        => Dist prob a -- ^
        -> a
        -> prob
(?=) p a = just a ??? p

liftC2 op p_ygx a b = op a (p_ygx b)

-- | Create a joint distribution of two independent discrete distributions
(<*>) :: (Num prob, Ord a, Ord b)
         => Dist prob a -- ^
         -> Dist prob b
         -> Dist prob (a,b)
(<*>) = liftM2 (,)

-- | Create a joint distribution of a distribution
-- and a conditional distribution
(-|-) :: (Num prob, Ord a, Ord b)
         => (a -> Dist prob b)
         -> Dist prob a -- ^
         -> Dist prob (b,a)
(-|-) p q = q >>= (\a -> liftM (flip (,) a) $! p a)

-- | Create an conditional joint distribution
(<\/>) :: (Num prob, Ord a, Ord b, Ord c)
          => (a -> Dist prob b)
          -> (a -> Dist prob c)
          -> (a -> Dist prob (b,c))
(<\/>) p q a = (p a) <*> (q a)


-- | Parallel conditional distributions
(<=>) :: (Num prob, Ord a, Ord b, Ord c, Ord d)
          => (a -> Dist prob b)
          -> (c -> Dist prob d)
          -> ((a,c) -> Dist prob (b,d))
(<=>) p q (a,c) = (p a) <*> (q c)

-- | Concatenate conditional distributions
(-<-) :: (Num prob, Ord a, Ord b, Ord c)
         => (b -> Dist prob c)
         -> (a -> Dist prob b)
         -> (a -> Dist prob (c,b))
(-<-) p q = liftC2 (-|-) q p

-- | Create a joint distribution of two independent discrete distributions
(<:>) :: (Num prob, Ord a) =>
          Dist prob a -- ^
          -> Dist prob [a]
          -> Dist prob [a]
(<:>) = liftM2 (:)

(<**>) :: (Num prob, Ord a) =>
          Dist prob [a] -- ^
          -> Dist prob [a]
          -> Dist prob [a]
(<**>) = liftM2 (++)

-- | Create a joint distribution of a distribution and a conditional distribution
(-|*-) :: (Num prob, Ord a)
         => ([a] -> Dist prob [a])
         -> Dist prob [a] -- ^
         -> Dist prob [a]
(-|*-) p q = q >>= (\a -> liftM (flip (++) a) $! p a)

-- | Create an conditional joint distribution
(<\/*>) :: (Num prob, Ord a)
          => ([a] -> Dist prob [a])
          -> ([a] -> Dist prob [a])
          -> ([a] -> Dist prob [a])
(<\/*>) p q a = (p a) <**> (q a)

-- | Create a parallel joint distribution
(<=*>) :: (Num prob, Ord a)
          => (a -> Dist prob [a]) -- ^ Singleton argument
          -> ([a] -> Dist prob [a])
          -> ([a] -> Dist prob [a])
(<=*>) p q (a:as) = (p a) <**> (q as)

(<=:>) :: (Num prob, Ord a)
          => (a -> Dist prob a) -- ^ Singleton argument
          -> ([a] -> Dist prob [a])
          -> ([a] -> Dist prob [a])
(<=:>) p q (a:as) = (p a) <:> (q as)

-- | Concatenate conditional distributions
(-<*-) :: (Num prob, Ord a)
          => ([a] -> Dist prob [a])
          -> ([a] -> Dist prob [a])
          -> ([a] -> Dist prob [a])
(-<*-) p q = liftC2 (-|*-) q p

-- | Reverse the conditional distribution
(-\-) :: (Fractional prob, Eq prob, Ord a, Ord b)
         => (a -> Dist prob b)
         -> Dist prob a -- ^
         -> (b -> Dist prob a)
(-\-) p q =
  let r = p -:- q
  in (\b -> fromFreqs' [(x,v * (p x ?= b) / (r ?= b))
                       | (x,v) <- decons q, v /= 0 && (p x ?= b) /=0])

-- | Concatenate conditional distributions and marginalize over b
(-.-) :: (Fractional prob, Ord a, Ord b, Ord c)
         => (b -> Dist prob c) -- ^
         -> (a -> Dist prob b)
         -> (a -> Dist prob c)
(p -.- q) a = fromFreqs' [(v*w) `seq` (c, v*w)
                         | (b,v) <- decons (q a), (c,w) <- decons (p b) ]

(-:-) :: (Num prob, Ord a, Ord b)
         =>  (a -> Dist prob b) -- ^
         -> Dist prob a
         -> Dist prob b
(-:-) p q = margin fst (p -|- q)

-- | A normalized combination of two distribution
combine :: (Fractional prob, Ord a)
           => (prob -> prob -> prob) -- ^
           -> Dist prob a
           -> Dist prob a
           -> Dist prob a
combine f p q = fromFreqs' [(x,f w (q ?= x)) | (x,w) <- decons p]

-- | Just strict verion of Numeric.Probability.Distribution.map
margin :: (Num prob, Ord a, Ord b)
          => (a->b) -- ^
          -> Dist prob a
          -> Dist prob b
margin f p = normS $ fmap f p

-- | Lift to array indices
arrayIndex :: (Num prob, Ord a)
              => Dist prob a -- ^
              -> Dist prob [a]
arrayIndex = fmap (: [])

-- | Convex combination
convex :: (Num a)
          => a -- ^ t
          -> a -- ^ p
          -> a -- ^ q
          -> a -- ^ t*p + (1-t)*q
convex t p q = t*p + (1-t)*q

readDistribution :: (Read prob, Read a, Fractional prob, Ord a)
                    => String
                    -> Dist prob a
readDistribution str =
  let lns = lines str
      header = head lns
      dist = tail lns
  in case header of
    "C" -> error "Trying to read in conditional distribution"
    "D" -> let vars = read $ head dist :: Int
           in (fromFreqs' . List.map (readNewFormat vars)) (tail $ tail dist)
    _ -> (fromFreqs' . List.map readOldFormat) dist
  where readOldFormat l =
          let (w:ws) = words l
          in (read w, read $ head ws)
        readNewFormat v l =
          let w = words l
          in (read $ toListStr (take v w),read $ last w)

readCondDistribution :: (Read prob, Read a, Fractional prob, Ord a)
                        => String
                        -> (a -> Dist prob a)
readCondDistribution str =
  let lns = List.filter (not . List.null) $ lines str
      header = head lns
      dist = tail lns
  in case header of
    "D" -> error "Trying to read in distribution"
    "C" -> let vars = words $ head dist
               varsY = read $ head vars :: Int
               varsX = read $ last vars :: Int
               jd = Map.fromListWith (++) . List.map
                  (\l -> let (y,r) = splitAt varsY (words l)
                             (x,r2) = splitAt varsX (tail r)
                             p = head $ tail r2
                         in (read $ toListStr x, [(read $ toListStr y, read p)]))
                   $ tail $ tail dist
           in (\a -> fromFreqs' (fromMaybe [] $ Map.lookup a jd))
    _ -> let jd = Map.fromListWith (++) . List.map
                  (\l -> let (w:ws) = words l
                         in (read w, [(read $ head ws, read $ ws !! 1)])) .
                  tail . lines $ str
         in (\a -> fromFreqs' (fromMaybe [] $ Map.lookup a jd))

toListStr s = "[" ++ List.intercalate "," s ++ "]"


