function colormaps=customColormaps()

% Distinguishable plot colors
customColors(1,:)=[0,0,0];
customColors(end+1,:)=[230, 25, 75];
customColors(end+1,:)=[60, 180, 75];
customColors(end+1,:)=[0, 130, 200];
customColors(end+1,:)=[245, 130, 48];
customColors(end+1,:)=[145, 30, 180];
customColors(end+1,:)=[70, 240, 240];
customColors(end+1,:)=[210, 245, 60];
customColors(end+1,:)=[250, 190, 190];
customColors(end+1,:)=[0, 128, 128];
customColors(end+1,:)=[230, 190, 255];
customColors(end+1,:)=[170, 110, 40];
customColors(end+1,:)=[128, 0, 0];
customColors(end+1,:)=[170, 255, 195];
customColors(end+1,:)=[128, 128, 0];
customColors(end+1,:)=[255, 215, 180];
customColors(end+1,:)=[0, 0, 128];
customColors(end+1,:)=[128, 128, 128];
customColors(end+1,:)=[255, 255, 25];
customColors=customColors./255;
colormaps.plot=customColors;

load('cmaps.mat');
colormaps.imgColor=cmaps.imgColor;
colormaps.twoType=cmaps.twoType;

end