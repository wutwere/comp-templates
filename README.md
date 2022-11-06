# precompiling bits/stdc++.h
find bits/stdc++ with `g++ -H file.cpp` then cd into directory
and then `g++ stdc++.h -O2 -DLOCAL` or whatever flags
# font
[install font here](https://github.com/powerline/fonts/blob/master/UbuntuMono/Ubuntu%20Mono%20derivative%20Powerline.ttf)
right click topbar and go to properties to set
# vim-plug install
first make sure you have [git](https://git-scm.com/download/win)
```powershell
iwr -useb https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim |`
    ni "$(@($env:XDG_DATA_HOME, $env:LOCALAPPDATA)[$null -eq $env:XDG_DATA_HOME])/nvim-data/site/autoload/plug.vim" -Force
```
# installing Neovim
[get chocolatey](https://chocolatey.org/install#individual)
then do `choco install neovim -y`
now nvim.exe should be in `C:\tools\neovim\nvim-win64\bin`
