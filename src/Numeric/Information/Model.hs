{-#LANGUAGE FlexibleInstances #-}
{- |
Module      :  Numeric.Information.Model
Description :  Build graphical models
Copyright   :  (c) Malte Harder
License     :  MIT

Maintainer  :  malte.harder@gmail.com
Stability   :  experimental
Portability :  portable

-}

module Numeric.Information.Model
       (
         -- * Data Types
         ModelNode (..)
       , EdgeBinding
       , NodeDist (..)
       , ModelEdge
       , Model
       , NodeLabels
       , ModelDistribution
       , ModelCondDistribution
       , ModelError (..)
         -- * Model Construction
       , mkModel
       , reduceModel
         -- * Marginalization & Conditionalization
       , marginalize
       , conditionalize
         -- * Distribution Accessors
       , extract
       , extractC
         -- * Utility Functions
       , sourceNodes
       , condNodes
       , sinkNodes
       , autoLabel
       )
       where

import Numeric.Information.Distribution
import Numeric.Information.IT

import Data.Graph.Inductive
import Data.List
import Data.Either
import Data.Maybe

import Control.Monad
import Control.Monad.Instances

import qualified Numeric.Probability.Distribution as D

data ModelNode prob a = ModelNode
                        { label :: String
                        , dist :: (NodeDist prob a)
                        } deriving (Show)

data NodeDist prob a = SrcDist (Dist prob a)
                     | CondDist ([a] -> Dist prob a)

instance (Show a, Show prob, Num prob, Ord prob, Ord a) => Show (NodeDist prob a) where
  show (SrcDist d) = show d
  show (CondDist a) = "<Conditional Distribution>"

instance (Show a, Num prob, Ord prob, Ord a) => Show (a -> Dist prob a) where
  show _ = "<Conditional Distribution>"

type EdgeBinding = Int
type ModelEdge = (String, String, EdgeBinding)

type Model prob a = Gr (ModelNode prob a) EdgeBinding

type NodeLabels = [(String, Int)]

type ModelDistribution prob a = (NodeLabels, Dist prob [a])
type ModelCondDistribution prob a = (NodeLabels , NodeLabels,
                                     [a] -> Dist prob [a])

data ModelError = StdE
                | EdgeDescriptionE
                | CycleFoundE
                | EdgeBindingE
                | DistE
                  deriving (Show, Eq)

mkModel :: (Num prob, Ord a)
           => [ModelNode prob a]
           -> [ModelEdge]
           -> Either ModelError (Model prob a)
mkModel uNodes uEdges =
  let nodes = autoLabel uNodes
      edges = catMaybes . map (\(dom,cod,bind) ->
                                maybeEdge (firstNode dom nodes,
                                           firstNode cod nodes,
                                           bind)) $ uEdges
  in if (length edges /= (length uEdges)) then
       Left EdgeDescriptionE
     else
       checkModel $ mkGraph nodes edges



reduceModel :: (Fractional prob, Ord a)
               => Model prob a
               -> ModelDistribution prob a
reduceModel m =
  let src = sort $ sourceNodes m
      sink = sort $ sinkNodes m
      srcD = map (srcDist . (nodeDist m)) src
      jointSrcD = foldr1 (<**>) $ map arrayIndex srcD
      (offset, nodeDesc,jointCondD) = graphFold m sink
  in (nodeDesc ++ (map (\i -> (label $ fromJust $ lab m $ src !! i,
                               i + offset)) [0.. (length src)-1]) ,
      jointCondD -|*- jointSrcD)

marginalize :: (Fractional prob, Ord a)
               => [String]
               -> ModelDistribution prob a
               -> ModelDistribution prob a
marginalize vars (vmap, d) =
  (varAutoIndex vars,
   margin (filterIndex vars vmap) d)
  where filterIndex vars vmap vs =
          [ vs !! (fromJust $ lookup v vmap)  | v <- vars]

conditionalize :: (Fractional prob, Eq prob, Ord a)
                  => [String]
                  -> [String]
                  -> ModelDistribution prob a
                  -> ModelCondDistribution prob a
conditionalize vars kvars md =
  let (_, jointD) = marginalize (vars ++ kvars) md
      (_, kMarginD) = marginalize kvars md
      (_, vMarginD) = marginalize vars md
      vs = length vars
  in (varAutoIndex vars, varAutoIndex kvars,
   \k -> fromFreqs' [(x, if (kMarginD ?= k) == 0 then
                           0
                         else (jointD ?= (x ++ k)) / (kMarginD ?= k))
                    | (x,v) <- D.decons vMarginD ] )

extract :: ModelDistribution prob a -> Dist prob [a]
extract = snd

extractC :: ModelCondDistribution prob a -> ([a] -> Dist prob [a])
extractC (_,_,d) = d

graphFold :: (Fractional prob, Ord a)
             => Model prob a
             -> [Node]
             -> (Int, [(String,Int)], [a] -> Dist prob [a])
graphFold m ns =
  let nodesD = map (\n -> if (indeg m n) == 0 then
                             delta
                          else
                            (\a -> arrayIndex $ (condDist $ (nodeDist m) n) a)
                   ) ns
      conNodes = filter (\n -> indeg m n /= 0) ns
      paNodes = sort $ nub $ concatMap (\n -> if (indeg m n) == 0 then
                                                [n]
                                              else pre m n) $ ns
      numCParents = length $ filter (\n -> indeg m n /= 0) paNodes
      edgesL = map (\n -> if (indeg m n) == 0 then
                               [(n,n,0)]
                          else inn m n) ns
      curD = (\a ->
               foldr1 (<**>) $
               map (\(e,d) -> d (permute (mapEdges paNodes e) a)) $
               zip edgesL nodesD)
      curL = (map (\i -> (label $ fromJust $ lab m $ conNodes !! i,
                               fromJust $ elemIndex (conNodes !! i) ns))
              [0.. (length conNodes)-1])
  in if numCParents > 0 then
        let (offset, preLabels, preDist) = (graphFold m paNodes)
            in (offset + (length ns),
                shiftIndices (length ns) preLabels ++ curL,
                curD -<*- preDist)
     else (length ns,curL,curD)
  where mapEdges pa e =
          map (\(dom,cod,bind) -> (bind,fromJust $ elemIndex dom pa)) e
        shiftIndices k idcs =
          map (\(l,i) -> (l,i+k)) idcs

sourceNodes m = (filter (\n -> indeg m n == 0) . nodes) m
condNodes m = (filter (\n -> indeg m n /= 0) . nodes) m
sinkNodes m = (filter (\n -> outdeg m n == 0) . nodes) m

nodeDist m n =
  let nodes = labNodes m
  in dist $ fromJust $ lookup n nodes

srcDist (SrcDist d) = d
srcDist (CondDist _) = error "No source distribution in node"

condDist (SrcDist _) =  error "No conditional distribution in node"
condDist (CondDist d) = d

checkModel :: (Num prob, Ord a)
              => Model prob a
              -> Either ModelError (Model prob a)
checkModel m =
  do when ((length $ scc m) /= (length $ nodes m))
       (Left CycleFoundE)
     when ((length $ filter (checkBinding m) (nodes m)) > 0 )
       (Left EdgeBindingE)
     when ((length $ filter (isSource m ) (condNodes m)) > 0 )
       (Left DistE)
     when ((length $ filter (not . (isSource m)) (sourceNodes m)) > 0 )
       (Left DistE)
     return m

checkBinding :: (Num prob, Ord a)
                => Model prob a
                -> Node
                -> Bool
checkBinding m n =
  let edges = inn m n
  in (indeg m n) /=
     (length . nub $
      filter (\l -> l >= 0 && l < (indeg m n)) $
      map (\(_,_,l) -> l) edges)

isSource :: (Num prob, Ord a)
            => Model prob a
            -> Node
            -> Bool
isSource m n =
  let nodes = labNodes m
      d = dist $ fromJust $ lookup n nodes
  in case d of
    (SrcDist _) -> True
    (CondDist _) -> False

firstNode :: String -> [(Node, ModelNode prob a)] -> Maybe Node
firstNode lbl = (liftM fst) . find (\(n,mn) -> label mn == lbl)

maybeEdge :: (Maybe Node, Maybe Node, EdgeBinding) ->
             Maybe (Node, Node, EdgeBinding)
maybeEdge (dom,cod, bind) =
  do jdom <- dom
     jcod <- cod
     return (jdom,jcod, bind)

varAutoIndex vars = [(vars !! i,i) | i <- [0 .. (length vars)-1]]

permute :: [(Int,Int)] -> [a] -> [a]
permute pi xs = map (\i -> xs !! (fromJust $ lookup i pi))  [0..(length pi) - 1]

autoLabel :: [a] -> [(Node,a)]
autoLabel xs = zip ([0..(length xs)-1]) xs


