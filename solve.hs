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
--autoDeriv :: Num a => (D a -> D a) -> a -> a
--autoDeriv func x = deriv
  --where (D _ deriv) = func (D x 1)

{-
-- Implement Automatic Differentiation
-- From: http://conal.net/blog/posts/what-is-automatic-differentiation-and-why-does-it-work
data D a = D a a deriving (Eq,Show)

instance Num a => Num (D a) where
  D x x' + D y y' = D (x+y) (x'+y')
  D x x' * D y y' = D (x*y) (y'*x + x'*y)
  fromInteger x   = D (fromInteger x) 0
  negate (D x x') = D (negate x) (negate x')
  signum (D x _ ) = D (signum x) 0
  abs    (D x x') = D (abs x) (x' * signum x)

instance Fractional x => Fractional (D x) where
  fromRational x  = D (fromRational x) 0
  recip  (D x x') = D (recip x) (- x' / sqr x)

sqr :: Num a => a -> a
sqr x = x * x

instance Floating x => Floating (D x) where
  pi              = D pi 0
  exp    (D x x') = D (exp    x) (x' * exp x)
  log    (D x x') = D (log    x) (x' / x)
  sqrt   (D x x') = D (sqrt   x) (x' / (2 * sqrt x))
  sin    (D x x') = D (sin    x) (x' * cos x)
  cos    (D x x') = D (cos    x) (x' * (- sin x))
  asin   (D x x') = D (asin   x) (x' / sqrt (1 - sqr x))
  acos   (D x x') = D (acos   x) (x' / (-  sqrt (1 - sqr x)))
  atan   (D x x') = D (atan   x) (x' / (sqr x + 1))
  sinh   (D x x') = D (sinh   x) (x' * cosh x)
  cosh   (D x x') = D (cosh   x) (x' * sinh x)
  asinh  (D x x') = D (asinh  x) (x' / sqrt (sqr x + 1))
  atanh  (D x x') = D (atanh  x) (x' / (1 - sqr x))
  acosh  (D x x') = D (acosh  x) (x' / sqrt ((x - 1)*(x + 1)))
-}

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
  acosh   = chain acosh (\x -> - recip (sqrt (sqr x - 1)))
  atanh   = chain atanh (\x -> recip (1 - sqr x))
