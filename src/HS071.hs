module HS071
( getNumVars
, getNumCombinations
, getNumBinaryVariables
, getBinaryIndices
, createBinaryInputCombinations
, solve
)
where

import Foreign
import Foreign.Ptr
import Foreign.C.Types
import Foreign.Marshal.Array
import System.IO.Unsafe
import Data.List.Split

foreign import ccall unsafe "hs071_c.h solve"
     c_solve :: Ptr CDouble -> IO (CDouble)

foreign import ccall unsafe "hs071_c.h create_binary_input_combinations"
     c_createBinaryInputCombinations :: IO (Ptr CDouble)

foreign import ccall unsafe "hs071_c.h get_num_combinations"
     c_getNumCombinations :: CSize

foreign import ccall unsafe "hs071_c.h get_num_binary_variables"
     c_getNumBinaryVariables :: CSize

foreign import ccall unsafe "hs071_c.h get_num_vars"
     c_getNumVars :: CSize

foreign import ccall unsafe "hs071_c.h get_binary_indices"
     c_getBinaryIndices :: IO (Ptr CInt)


getNumVars :: Int
getNumVars = fromIntegral c_getNumVars


getNumCombinations :: Int
getNumCombinations = fromIntegral c_getNumCombinations


getNumBinaryVariables :: Int
getNumBinaryVariables = fromIntegral c_getNumBinaryVariables


getBinaryIndices' :: IO ([CInt])
getBinaryIndices' = do
            cdata <- c_getBinaryIndices
            let numBinaryVariables = getNumBinaryVariables
            peekArray numBinaryVariables cdata


getBinaryIndices :: [Int]
getBinaryIndices =
        let coutput = System.IO.Unsafe.unsafePerformIO getBinaryIndices'
            output  = map fromIntegral coutput
        in  output


createBinaryInputCombinations' :: IO ([CDouble])
createBinaryInputCombinations' = do
                        cdata <- c_createBinaryInputCombinations
                        let numElements = getNumBinaryVariables * getNumCombinations
                        peekArray numElements cdata


createBinaryInputCombinations :: [Double]
createBinaryInputCombinations =
            let coutput = System.IO.Unsafe.unsafePerformIO createBinaryInputCombinations'
                output = map realToFrac coutput
            in  output


-- not really infinity
infinity :: Double
infinity = 2e19 :: Double


solve' :: [Double] -> IO (Maybe Double)
solve' input = do
            let cinput = map realToFrac input
            optimalValue <- withArray cinput (\x -> do
                                    c_solve x)
            if  optimalValue > (realToFrac infinity) then
                return Nothing
                else return $ Just (realToFrac optimalValue)

solve :: [Double] -> Maybe Double
solve input = System.IO.Unsafe.unsafePerformIO $ solve' input

