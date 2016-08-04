%% Tutorial 2: Introducing the Link Transmission Model (LTM)

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
% This tutorial introduces a basic implementation of the link transmission 
% model. This module propagates traffic over a network guided by predefined
% turning fractions (or splitting rates) at nodes. Congestion and spillback
% is modelled with three different queuing principles. The standard queuing
% approach (most realistic for traffic modelling) is based on first-order 
% kinematic wave theory and a triangular fundamental diagram. Implemented 
% alternatives are vertical queueing and horizontal queueing. 
%

%add these folders to the search path
addpath('Dynamic Traffic Assignment','Visualization Tools','Network Data')
%clear the work space
clear
%clear the command window
clc
%close all windows
close all

display('<<<Introducing the Link Transmission Model (LTM)>>>')

%% Loading the data
% The network represents a two-lane highway with vehicles moving from right 
% to left. There are two on-ramps feeding additional traffic into the 
% system. The demand pattern is chosen such that the most downstream
% merge forms a temporary bottleneck.
%

% Network and demand data
load net1.mat

% Plot the network
plotNetwork(nodes,links,true,[]);

%% Setup the simulation
% Before the simulation can be run the time interval has to be set and the
% total number of time steps has to be defined. These are used to transform
% the different origin-destination (OD-) matrices into a 3D-matrix. The OD
% matrices represent stationary demand rates inbetween the time tics
% defined in the timeSeries variable. This makes it easy to adjust the time
% interval (dt) of the simulation and update the input ODmatrix by 
% rerunning this section. The turning fractions have to be defined for each 
% time interval. They are stored in a cell array where each cell represents
% the turns in matrix with size [#incoming links, #outgoing links]. All 
% turning fractions are equal to one because no diverges are present in 
% this network.
%

%setup the time interval and total number of time steps
dt = 0.002;
totT = 2/dt;

%build the full ODmatrix
[ODmatrix,origins,destinations] = buildODmatrix(ODmatrices,timeSeries,dt,totT);

%initilize the Turning Fractions
TF = num2cell(ones(size(nodes.id,1),totT));
for t=1:totT
    TF{10,t} = ones(2,1);
    TF{25,t} = ones(2,1);
end
clear t;

%% Visualize the demand in the network
% The demand between every origin-destination combination is plotted for
% each time interval of the simulation.
%

figure;
plot(dt*[0:totT-1],reshape(ODmatrix(1,1,:),1,[]),'d-b',dt*[0:totT-1],reshape(ODmatrix(2,1,:),1,[]),'.-r',dt*[0:totT-1],reshape(ODmatrix(3,1,:),1,[]),'x-g')
legend('OD 1-28','OD 31-28','OD 34-28')
xlabel('time (h)')
ylabel('flow (veh/h)')

%% Compute the single-commodity Dynamic Network Loading
% The link transmission model propagates the traffic in the simulation. It
% is composed of link model and node model which are both consistent with 
% first-order traffic flow theory with a triangular fundamental diagram.
% For each link it is therefor required to define free-flow speed, capacity
% and jam density. In hypercritical traffic states this results in delays
% that are function of the density of vehicles. Alternatively it is also
% possible to simulate vertical queues that do not propagate upstream or
% horizontal queues that characterized by a single density for all 
% hypercritical states. In this implementation the queue principle is set 
% identical for all links in the network by suppling a 2-character string
% to the LTM_SC function:
%
% * PQ: physical queue (standard)
% * VQ: vertical queue
% * HQ: horizontal queue
%
% Within the LTM_SC function the *sending flow* (representing the demand
% towards a node) and the *receiving flow* (representing the maximum 
% number of vehicles alowed to propagate into a link) are evaluated for 
% every node in each time interval. The nested functions that compute them 
% are found at the end of the script. The node model governs the transfer 
% of flow over a node and is included as a seperate function.
%

display('Running LTM')

%link model
mode = 'PQ'; 
%PQ: physical queue (standard)
%VQ: vertical queue
%HQ: horizontal queue

%run LTM
tic
[cvn_up,cvn_down] = LTM_SC(nodes,links,origins,destinations,ODmatrix,dt,totT,TF,mode);
toc

%% Visualize the resulting densities and flows using XT diagrams
% The result of the link transmission model is expressed in cumulative 
% vehicle numbers (CVN) for every links upstream and downstream end over 
% the time domain. Flows and densities are computed in post processing 
% phase as time and space derivates of these CVN functions. Density and 
% flow are depicted in space-time (or XT) diagrams of the main road.
%

%compute the simulated densities & flows
[simDensity] = cvn2dens(cvn_up,cvn_down,totT,links);
[simFlows_down] = cvn2flows(sum(cvn_down,3),dt);

%Main road
plotXT(links,1:27,simDensity,dt,totT);
title('XT-graph of densities on the main road: LTM','FontSize',14,'fontweight','b')
plotXT(links,1:27,simFlows_down,dt,totT-1);
title('XT-graph of downstream flows on the main road: LTM','FontSize',14,'fontweight','b')

%% Transform CVN values to travel times
% The upstream and dowsntream CVN functions are transformed into travel 
% times for every link in the network. This requires computing the time 
% window between the occurence of a specific vehicle number on upstream and
% downstream end of the link.
%

%calculate the simulated travel times
[simTT] = cvn2tt(cvn_up,cvn_down,dt,totT,links);

%visualize the travel time along the main route (from split to merge)
[~,~,~,tt]=plotTT(links,1:27,simTT,dt,totT);

%% Make an animation of the result
% The variation of densities is animated in the network by considering only 
% every 10th simulation interval.
%

% fRate = 20; %set frame rate
fRate = inf; %allows the for manual control using space bar
animateSimulation(nodes,links,simDensity(:,1:10:end),dt*[0:10:totT],fRate); %only shows every 10th frame

%% Closing notes
%
% * The network is deliberatly chosen to be very detailed. This is done for
% visualizing the results more accuratly. In tutorial 4 a coarser network
% is presented with the same layout.
% * To compare the queuing principles it is required to rerun the script
% while the mode parameter is adjusted.
% 