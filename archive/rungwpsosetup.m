function rungwpsosetup(paramsFile,outdataFilePrfx)
rungwpso(paramsFile,outdataFilePrfx) %pwelch 
paramsFileshps = [paramsFile,'shps'];
rungwpso(paramsFileshps,[outdataFilePrfx,'shps_']) %shapes estimate
end