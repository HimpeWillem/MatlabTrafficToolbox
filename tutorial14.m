%% Tutorial 14: Advanced deterministic equilibrium routing with warm start

%% Disclaimer
% This file is part of the matlab package for dynamic traffic assignments 
% developed by the KULeuven. 
%
% Copyright (C) 2016  Himpe Willem, Leuven, Belgium
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% More information at: http://www.mech.kuleuven.be/en/cib/traffic/downloads
% or contact: willem.himpe {@} kuleuven.be

%% Introduction
% In this tutorial an advanced algorithm is implemented to find 
% deterministic user equilibrium in a dynamic traffic assignment. The 
% solution is found by means of reduced gradient projection which is 
% applied to all nodes simultaneously or sequentially.
%

%add these folders to the search path
addpath('Dynamic Traffic Assignment','Visualization Tools','Network Data')
javaclasspath('Dynamic Traffic Assignment');
%clear the work space
clear
%clear the command window
close all

display('<<<Advanced deterministic equilibrium routing with warm start>>>')
