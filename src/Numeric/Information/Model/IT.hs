{- |
Module      :  Numeric.Information.Model.IT
Description :  Information quantities on models
Copyright   :  (c) Malte Harder
License     :  MIT

Maintainer  :  malte.harder@gmail.com
Stability   :  experimental
Portability :  portable

-}

module Numeric.Information.Model.IT
       (
         -- * Entropy & Mutual Information
         nodeEntropy
       , nodeCondEntropy
       , nodeMutualInfo
       , nodeCondMutualInfo
         -- * Specific Information
       , nodeSpecificInformation
       , nodeMinimalInformation
         -- * Partial Information Decomposition
         -- ** Data Types
       , PINode (..)
       , PILattice
         -- ** Functions
       , piDecomp
         -- ** Helper Functions
       , piLattice
       , piSet
       , (-<<)
       , (-<=)
       , subset
         -- * Information Flow
       , nodeInformationFlow
       ) where

import Numeric.Information.Model

import Numeric.Information.Distribution
import Numeric.Information.IT

import Data.Graph.Inductive
import Data.List

--import Statistics.Math hiding(log2)
import Numeric.SpecFunctions hiding(log2)


nodeEntropy :: (Floating prob, Ord prob, Ord a)
               => [String]
               -> ModelDistribution prob a
               -> prob
nodeEntropy v m = (entropy . extract) (marginalize v m)

nodeCondEntropy :: (Floating prob, Ord prob, Ord a)
                   => [String]
                   -> [String]
                   -> ModelDistribution prob a
                   -> prob
nodeCondEntropy yv xv m =
  let p_ygx = extractC $ conditionalize yv xv m
      p_x = extract $ marginalize xv m
  in condEntropy p_ygx p_x

nodeMutualInfo :: (Floating prob, Ord prob, Ord a)
                   => [String]
                   -> [String]
                   -> ModelDistribution prob a
                   -> prob
nodeMutualInfo yv xv m =
  let p_ygx = extractC $ conditionalize yv xv m
      p_x = extract $ marginalize xv m
  in mutualInfo p_ygx p_x

nodeCondMutualInfo :: (Floating prob, Ord prob, Ord a)
                      => [String]
                      -> [String]
                      -> [String]
                      -> ModelDistribution prob a
                      -> prob
nodeCondMutualInfo yv xv zv m =
  let p_xgz = extractC $ conditionalize xv zv m
      p_xgyz = \(y,z) -> (extractC $ conditionalize xv (yv++zv) m) $ y ++ z
      p_z = extract $ marginalize zv m
      p_ygz = extractC $ conditionalize yv zv m
  in condMutualInfo p_xgz p_xgyz p_z (p_ygz -|- p_z)

nodeSpecificInformation :: (Floating prob, Ord prob, Ord a)
                           => [String]
                           -> [a]
                           -> [String]
                           -> ModelDistribution prob a
                           -> prob
nodeSpecificInformation sv s av m =
  let p_s = (extract $ marginalize sv m) ?= s
      p_ags = (extractC $ conditionalize av sv m) s
      p_sga = (extractC $ conditionalize sv av m)
  in expected p_ags (\a-> log2 ((p_sga a) ?= s / p_s))

nodeMinimalInformation :: (Floating prob, Ord prob, Ord a)
                           => [String]
                           -> [[String]]
                           -> ModelDistribution prob a
                           -> prob
nodeMinimalInformation  sv avs m =
  let p_s = (extract $ marginalize sv m)
  in expected p_s (\s -> minimum $
                         map (\av -> nodeSpecificInformation sv s av m) avs)

nodeMaxMinimalInformation :: (Floating prob, Ord prob, Ord a)
                           => [String]
                           -> [[[String]]]
                           -> ModelDistribution prob a
                           -> prob
nodeMaxMinimalInformation sv avss m =
  let p_s = (extract $ marginalize sv m)
  in expected p_s (\s -> maximum' $
                         map (\avs -> minimum $
                                      map (\av -> nodeSpecificInformation sv s av m)
                                      avs) avss)


nodeInformationFlow :: (Floating prob, Ord prob, Ord a)
                       => [String]
                       -> [String]
                       -> Model prob a
                       -> prob
nodeInformationFlow as bs m = undefined


maximum' [] = 0
maximum' xs = maximum xs

data PINode prob a = PINode { partial :: [[a]]
                            , value :: prob
                            } deriving (Eq)

instance (Show a, Show prob) => Show (PINode prob a) where
  show n = (intercalate "," $
            map (\xs -> (intercalate " & " (map (trim . show . show) xs))) $ partial n)
           ++ " = "++ (show $ value n)

trim = reverse . tail . reverse . tail

type PILattice prob a = Gr (PINode prob a) ()

piDecomp :: (Floating prob, Ord prob, Ord a)
            => [String]
            -> [String]
            -> ModelDistribution prob a
            -> PILattice prob String
piDecomp s as m =
  let g  = piLattice as
      g' = nmap minimalInformation g
  in gmap (subtractLower g') g'
  where minimalInformation piN =
          piN { value = (nodeMinimalInformation s (partial piN) m)}
        subtractLower g' (inA, n, piN, outA) =
          (inA,
           n,
           piN { value = (value piN) -
                         (nodeMaxMinimalInformation s
                          (map (partial . snd) $ filter
                          (\(n',piN') -> (n',n) `elem` (edges g') ) $ labNodes g')
                          m )} ,
           outA)

piLattice :: (Eq a, Fractional prob) => [a] -> PILattice prob a
piLattice modelNodes =
  let nodes = autoLabel $ map (\n -> PINode { partial = n, value = 0 })
              $ piSet modelNodes
      size = (length nodes)-1
      edges = [ (i,j,()) | i <- [0..size], j <- [0..size],
                i /= j && -- Reflexive reduction
                (partial $ snd $ nodes !! i) -<= (partial $ snd $ nodes !! j) ]
      compositionEdges = [ (i,k,()) | (i,j,()) <- edges, k <- [0..size],
                               (j,k,()) `elem` edges]
      reducedEdges = edges \\ compositionEdges
  in mkGraph nodes reducedEdges

piSet :: (Eq a) => [a] -> [[[a]]]
piSet r = [ a | a <- (powerset' . powerset') r,
          (null [ a | x <- a, y <- a,
                  x /= y &&
                  (length x <= length y) && (subset x y)])]

muSet :: (Eq a) => Int -> [a] -> [([[a]], Double)]
muSet k ls =
  let r = piSet ls
      p = powerset' ls
      n = length ls
  in map (\beta -> (beta, (fromIntegral n)/(fromIntegral k) * (esoohc n k) * (fromIntegral $ length $
                           filter (\a -> length a == k && beta -<= [a]) p) - 1 ) ) r


muSet' :: (Eq a) => [a] -> [([[a]], Double)]
muSet' ls =
  let r = piSet ls
      p = powerset' ls
      n = length ls
  in map (\beta ->
           (beta,
            (sum [ (esoohc n k)  * ( fromIntegral
                    $ (length $ filter
                   (\a -> length a == k && beta -<= [a]) p)) - (fromIntegral k)/(fromIntegral n)
                 | k <- [1 .. (n-1)] ] ) ) ) r

esoohc :: Int -> Int -> Double
esoohc n k = 1.0 / (n `choose` k)

(-<<) :: (Eq a) => [[a]] -> [[a]] -> Bool
(-<<) a b = (a -<= b) && (a /= b)

(-<=) :: (Eq a) => [[a]] -> [[a]] -> Bool
(-<=) a b = foldr q True b
  where q bel cur = cur && (not . null $ filter (\ael -> subset ael bel)  a )

subset :: (Eq a) => [a] -> [a] -> Bool
subset a b = foldr q True a
  where q el cur = cur && (elem el b)

powerset' :: [a] -> [[a]]
powerset' = tail . powerset

powerset       :: [a] -> [[a]]
powerset []     = [[]]
powerset (x:xs) = xss /\/ map (x:) xss
  where xss = powerset xs

(/\/)        :: [a] -> [a] -> [a]
[]     /\/ ys = ys
(x:xs) /\/ ys = x : (ys /\/ xs)