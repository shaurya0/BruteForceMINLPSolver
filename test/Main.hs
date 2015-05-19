module Main where

import HS071
import Data.List.Split
import Control.Parallel.Strategies
import Control.Parallel

bruteForceMinlpSolve :: Int -> Int -> [Double] -> [Maybe Double]
bruteForceMinlpSolve numBinaryVariables numCombinations inputCombinations =
                let inputCombinations = createBinaryInputCombinations
                    inputs = chunksOf numBinaryVariables inputCombinations
                    results = (parMap rseq) solve inputs
                in results

main = do
    let numBinaryVariables = getNumBinaryVariables
    let numCombinations    = getNumCombinations
    let inputCombinations  = createBinaryInputCombinations
    let results = bruteForceMinlpSolve numBinaryVariables numCombinations inputCombinations
    putStrLn $ show results
