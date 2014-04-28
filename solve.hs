--  Solve, a Haskell module that enables finding roots of equations
--  Copyright (C) 2014 Matthew Filmer
--
--  This program is free software: you can redistribute it and/or modify
--  it under the terms of the GNU General Public License as published by
--  the Free Software Foundation, either version 3 of the License, or
--  (at your option) any later version.
--
--  This program is distributed in the hope that it will be useful,
--  but WITHOUT ANY WARRANTY; without even the implied warranty of
--  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
--  GNU General Public License for more details.
--
--  You should have received a copy of the GNU General Public License
--  along with this program.  If not, see <http://www.gnu.org/licenses/>.

module Solve (
  search, 
  bisection, 
  ridder, 
  newton, 
  simpNewton, 
  numDeriv,
  Deriv (D)) where

-- search: Approximate the root of a function given a starting point
--         and an increment
search :: (Enum a, Num a, Num b, Ord b) => (a -> b) -> a -> a -> (a, a)
search func start inc =
  head [(a, b) | (a, b) <- zip vals (tail vals), func a * func b <= 0]
    where
      vals = [start, start + inc ..]

-- bisection: Calculate the root of a function given a starting range
bisection :: (Fractional a, Ord a, Num b, Ord b) => (a-> b) -> (a, a) -> a
bisection func (low, high)
  | abs (high - low) < abs (1e-14 * mid) = mid
  | otherwise = if func low * func mid <= 0
                  then bisection func (low, mid)
                  else bisection func (mid, high)
    where
      mid = (high + low) / 2

ridder :: (Floating a, Ord a) => (a -> a) -> (a, a) -> a
ridder func (low, high)
  | abs (high - low) < precision = mid
  | otherwise = if midV * rootV <= 0
                  then ridder func (mid, root)
                  else
                    if lowV * rootV < 0
                      then ridder func (low, root)
                      else ridder func (root, high)
    where
      mid = (high + low) / 2
      lowV = func low
      midV = func mid
      highV = func high
      rootV = func root
      root
        | lowV > highV = mid + (mid - low) * midV / sqrt (midV**2 - lowV*highV)
        | otherwise = mid - (mid - low) * midV / sqrt (midV**2 - lowV*highV)
      precision = max (abs mid * 1e-14) 1e-14

-- Find a root using Newton's Method, iterating until the root doesn't change
fullNewton :: (Eq a, Fractional a) => (a -> a) -> (a -> a) -> a -> a
fullNewton func deriv guess = approxs !! final
  where
    approxs = iterate (newton' func deriv) guess
    final = length $ takeWhile (\(x, y) -> x /= y) (zip approxs (tail approxs))

simpFullNewton :: (Fractional a, Ord a) => (a -> a) -> a -> a
simpFullNewton func guess = fullNewton func (numDeriv func) guess

newton :: (Fractional a) => (a -> a) -> (a -> a) -> a -> Int -> a
newton func deriv guess n = iterate (newton' func deriv) guess !! (n - 1)

simpNewton :: (Fractional a, Ord a) => (a -> a) -> a -> Int -> a
simpNewton func guess n = newton func (numDeriv func) guess n

newton' :: (Fractional a) => (a -> a) -> (a -> a) ->  a -> a
newton' func deriv guess = guess - ((func guess) / (deriv guess))

-- Numerically calculate the derivative of a function at a point
numDeriv :: (Fractional a, Ord a) => (a -> a) -> a -> a
numDeriv func x = (func (x + dx) - func (x - dx)) / (2 * dx)
  where
    dx = max 1e-10 (abs x * 1e-10)

-- Automatically differentiate (implementation below)

-- Derivative of an input value
dx x = D x (dConst 1)

-- First derivative
autoDeriv :: Num a => (Deriv a -> Deriv a) -> a -> a
autoDeriv = autoDerivN 1

-- Nth derivative
autoDerivN :: Num a => Int -> (Deriv a -> Deriv a) -> a -> a
autoDerivN n func x = getNthDeriv n (func (dx x))
  where
    getNthDeriv 0 (D x _) = x
    getNthDeriv n (D _ x') = getNthDeriv (n - 1) x'

-- List of derivative values
autoDerivAll :: Num a => (Deriv a -> Deriv a) -> a -> [a]
autoDerivAll func x = l (func (dx x)) where l (D a0 a') = a0:l a'

-- List of derivative functions
--allAutoDeriv :: Num a => (Deriv a -> Deriv a) -> [Deriv a -> Deriv a]
--allAutoDeriv func = l where l = l:allAutoDeriv l
--allAutoDeriv func = l:autoDeriv l where l = autoDeriv func
--allAutoDeriv func = ad:allAutoDeriv ad where ad = autoDeriv func

-- Implementation of Automatic Differentiation for an unlimited number of
-- derivitives
-- Based on http://conal.net/blog/posts/beautiful-differentiation

-- Datatype
data Deriv a = D a (Deriv a)

-- Derivative of a constant
dConst :: Num a => a -> Deriv a
dConst x0 = D x0 (dConst 0)

-- Implementation of chain rule
chain f f' u@(D u0 u') = D (f u0) (f' u * u')

-- Squaring function
sqr x = x*x

-- Basic operators
instance Num a =>  Num (Deriv a) where
  D x0 x' + D y0 y'           = D (x0 + y0) (x' + y')
  D x0 x' - D y0 y'           = D (x0 - y0) (x' - y')
  x@(D x0 x') * y@(D y0 y')   = D (x0 * y0) (x' * y + x * y')
  fromInteger                 = dConst . fromInteger
  negate                      = chain negate (\_ -> dConst (-1))
  abs                         = chain abs signum
  signum                      = chain signum (\_ -> dConst (0))

instance Fractional a => Fractional (Deriv a) where
  x@(D x0 x') / y@(D y0 y')   = D (x0 / y0) ((x' * y - x * y') / (sqr y))
  fromRational  = dConst . fromRational
  recip         = chain recip (negate . sqr . recip)

-- Functions
instance Floating a => Floating (Deriv a) where
  pi      = dConst pi
  exp     = chain exp exp
  log     = chain log recip
  sqrt    = chain sqrt (recip . (2*) . sqrt)
  sin     = chain sin cos
  cos     = chain cos (negate . sin)
  asin    = chain asin (\x -> recip (sqrt (1 - sqr x)))
  acos    = chain acos (\x -> - recip (sqrt (1 - sqr x)))
  atan    = chain atan (\x -> recip (sqr x + 1))
  sinh    = chain sinh cosh
  cosh    = chain cosh sinh
  asinh   = chain asinh (\x -> recip (sqrt (sqr x + 1)))
  acosh   = chain acosh (\x -> recip (sqrt (sqr x - 1)))
  atanh   = chain atanh (\x -> recip (1 - sqr x))

instance Show a => Show (Deriv a) where
  show (D a (D b _)) = "D " ++ show a ++ " (D " ++ show b ++ " (D ...))"
