%% Nhia & Nolen Belle

% Instructions: Follow through this code step by step, while also referring
% to the overall instructions and questions from the lab assignment sheet.

%% 1. Read in the monthly gridded CO2 data from the .csv file
% The data file is included in your repository as “LDEO_GriddedCO2_month_flux_2006c.csv”
% Your task is to write code to read this in to MATLAB
% Hint: you can again use the function “readtable”, and use your first data lab code as an example.
%<--
filename = 'LDEO_GriddedCO2_month_flux_2006c.csv'
CO2data = readtable(filename);
%% 2a. Create new 3-dimensional arrays to hold reshaped data
%Find each unique longitude, latitude, and month value that will define
%your 3-dimensional grid
longrid = unique(CO2data.LON); %finds all unique longitude values
latgrid = unique(CO2data.LAT); %<-- following the same approach, find all unique latitude values
monthgrid = unique(CO2data.MONTH); %<-- following the same approach, find all unique months

%Create empty 3-dimensional arrays of NaN values to hold your reshaped data
    %You can make these for any variables you want to extract - for this
    %lab you will need PCO2_SW (seawater pCO2) and SST (sea surface
    %temperature)
pCO2Array = NaN*zeros(length(longrid),length(latgrid), length(monthgrid));
sstArray = NaN*zeros(length(longrid),length(latgrid), length(monthgrid));



%% 2b. Pull out the seawater pCO2 (PCO2_SW) and sea surface temperature (SST)
%data and reshape it into your new 3-dimensional arrays

for i = 1:height(CO2data)
    PCO2_SW = CO2data.PCO2_SW(i);
    sst_SW = CO2data.SST(i);
    lat = CO2data.LAT(i);
    lon = CO2data.LON(i);
    month = CO2data.MONTH(i);
    indexLon = find(longrid == lon);
    indexLat = find(latgrid == lat);
    indexMonth = find(monthgrid == month);
    pCO2Array(indexLon,indexLat,indexMonth) = PCO2_SW;
    sstArray(indexLon,indexLat,indexMonth) = sst_SW;    
end   

%% 3a. Make a quick plot to check that your reshaped data looks reasonable
%Use the imagesc plotting function, which will show a different color for
%each grid cell in your map. Since you can't plot all months at once, you
%will have to pick one at a time to check - i.e. this example is just for
%January

imagesc(sstArray(:,:,1));
imagesc(pCO2Array(:,:,1));

%% 3b. Now pretty global maps of one month of each of SST and pCO2 data.
%I have provided example code for plotting January sea surface temperature
%(though you may need to make modifications based on differences in how you
%set up or named your variables above).

figure(1); clf
worldmap world
contourfm(latgrid, longrid, sstArray(:,:,1)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January Sea Surface Temperature (^oC)')


figure(2); clf
worldmap world
contourfm(latgrid, longrid, pCO2Array(:,:,1)','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('January pCO2 in Surface Water (uatm)')

%Check that you can make a similar type of global map for another month
%and/or for pCO2 using this approach. Check the documentation and see
%whether you can modify features of this map such as the contouring
%interval, color of the contour lines, labels, etc.

%<--

%% 4. Calculate and plot a global map of annual mean pCO2
%<-- 
meanpCO2 = mean(pCO2Array,3);
meanSST = mean(sstArray,3);
figure(3); clf
worldmap world
contourfm(latgrid, longrid, meanpCO2','linecolor','none');
colorbar
geoshow('landareas.shp','FaceColor','black')
title('Mean Annual pCO2 in Surface Water (uatm)')

%% 5. Calculate and plot a global map of the difference between the annual mean seawater and atmosphere pCO2
%<--


%% 6. Calculate relative roles of temperature and of biology/physics in controlling seasonal cycle
%<--
pCO2_BP = pCO2Array.^(0.0423*(repmat(meanSST,[1 1 12]))-sstArray);
pCO2_T = repmat(meanpCO2, [1 1 12]).^(0.0423*(sstArray - (repmat(meanSST,[1 1 12]))));

%% 7. Pull out and plot the seasonal cycle data from stations of interest
%Do for BATS, Station P, and Ross Sea (note that Ross Sea is along a
%section of 14 degrees longitude - I picked the middle point)
%at each location, we need the pCO2_BP, pCO2_T, and pCO2Array values over
%the season
BATSdata = NaN*zeros(12,3);
indexLonBATS = find(longrid == 312.5);
indexLatBATS = find(latgrid == 32);

RossData = NaN*zeros(12,3);
indexLonRoss = find(longrid == 177.5);
indexLatRoss = find(latgrid == -76);

Pdata = NaN*zeros(12,3);
indexLonP = find(longrid == 217.5);
indexLatP = find(latgrid == 48);
for i = 1:12
    BATSdata(i,1) = pCO2_BP(indexLonBATS, indexLatBATS,i);
    BATSdata(i,2) = pCO2_T(indexLonBATS, indexLatBATS,i);
    BATSdata(i,3) = pCO2Array(indexLonBATS, indexLatBATS,i);
    
    RossData(i,1) = pCO2_BP(indexLonRoss, indexLatRoss,i);
    RossData(i,2) = pCO2_T(indexLonRoss, indexLatRoss,i);
    RossData(i,3) = pCO2Array(indexLonRoss, indexLatRoss,i);
    
    Pdata(i,1) = pCO2_BP(indexLonP, indexLatP,i);
    Pdata(i,2) = pCO2_T(indexLonP, indexLatP,i);
    Pdata(i,3) = pCO2Array(indexLonP, indexLatP,i);
end 


% station = actual cordinates -> closest in the data
%BATS = 32.5, 296 -> 32,312.5
%Ross = -76.5, 176 -> -76,177.5
%Station P = 50, 216-> 48, 217.5
%% 8. Reproduce your own versions of the maps in figures 7-9 in Takahashi et al. 2002
% But please use better colormaps!!!
% Mark on thesese maps the locations of the three stations for which you plotted the
% seasonal cycle above

%<--
