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

module Solve (search, bisection, ridder, newton, simpNewton, numDeriv) where

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

simpNewton :: (Fractional a, Ord a) => (a -> a) -> a -> Int -> a
simpNewton func guess n = newton func (numDeriv func) guess n

newton :: (Fractional a) => (a -> a) -> (a -> a) -> a -> Int -> a
newton func deriv guess n = iterate (newton' func deriv) guess !! (n - 1)

newton' :: (Fractional a) => (a -> a) -> (a -> a) ->  a -> a
newton' func deriv guess = guess - ((func guess) / (deriv guess))

-- Numerically calculate the derivative of a function at a point
numDeriv :: (Fractional a, Ord a) => (a -> a) -> a -> a
numDeriv func x = (func (x + dx) - func (x - dx)) / (2 * dx)
  where
    dx = max 1e-14 (abs x * 1e-14)
