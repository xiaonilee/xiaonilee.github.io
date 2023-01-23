
# Installing rattle on MacOS 10.11 (or above)

# install gtk+ 2.24.32_3
system('brew install gtk+')


# macOS version:catalina 
local({
  if (Sys.info()[['sysname']] != 'Darwin') return()
  
  .Platform$pkgType = 'mac.binary.catalina'
  unlockBinding('.Platform', baseenv())
  assign('.Platform', .Platform, 'package:base')
  lockBinding('.Platform', baseenv())
  
  options(
    pkgType = 'both', install.packages.compile.from.source = 'always',
    repos = 'https://macos.rbind.io'
  )
})

install.packages(c('RGtk2', 'cairoDevice', 'rattle'))
