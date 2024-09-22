# Precompiling bits/stdc++.h
Create any .cpp file that includes <bits/stdc++.h>. Example:
```cpp
// a.cpp
#include <bits/stdc++.h>
using namespace std;
int main() {
    return 0;
}
```
Then, `g++ -H a.cpp` to find the directory of bits/stdc++.h, which should be near the top of the results.

Change into the directory you found and do `g++ -O2 -DLOCAL -std=c++17 stdc++.h` or your preferred flags

# Font
[Install font here](https://github.com/powerline/fonts/blob/master/UbuntuMono/Ubuntu%20Mono%20derivative%20Powerline.ttf), right click the top bar of vim/nvim, go to Properties->Font, set font

(9/21/2024 update: Cascadia Mono, which is on Windows by default I think, looks better)

# Installing vim-plug on Windows
First, make sure you have [Git](https://git-scm.com/download/win). Then, run this in Windows PowerShell:
```powershell
iwr -useb https://raw.githubusercontent.com/junegunn/vim-plug/master/plug.vim |`
    ni "$(@($env:XDG_DATA_HOME, $env:LOCALAPPDATA)[$null -eq $env:XDG_DATA_HOME])/nvim-data/site/autoload/plug.vim" -Force
```
# Installing Neovim on Windows
[First, install Chocolatey](https://chocolatey.org/install#individual).
Then, do `choco install neovim -y`, and now nvim.exe should be in `C:\tools\neovim\nvim-win64\bin`
