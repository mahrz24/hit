{- |
Module      :  Numeric.Information
Description :  Information theory framework for Haskell
Copyright   :  (c) Malte Harder
License     :  MIT

Maintainer  :  malte.harder@gmail.com
Stability   :  experimental
Portability :  portable


-}

module Numeric.Information (
    module Numeric.Information.IT
  , module Numeric.Information.Distribution
  , module Numeric.Information.Model
  , module Numeric.Information.Model.IT
  ) where

import Numeric.Information.Distribution
import Numeric.Information.IT
import Numeric.Information.Model
import Numeric.Information.Model.IT