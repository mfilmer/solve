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

module Solve where

search :: (Enum a, Num a, Num b, Ord b) => (a -> b) -> a -> a -> (a, a)
search func start inc =
  head [(a, b) | (a, b) <- zip vals (tail vals), (func a * func b) <= 0]
    where
      vals = [start, start + inc ..]