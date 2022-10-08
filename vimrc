call plug#begin()

Plug 'dense-analysis/ale'
let g:ale_linters = {'cpp': ['g++']}
let g:ale_cpp_cc_executable = 'g++'
let g:ale_cpp_cc_options = '-std=c++17 -Wall -Wextra -DLOCAL'

Plug 'vim-airline/vim-airline'
let g:airline_powerline_fonts = 1
let g:airline_theme = 'dracula'
let g:airline#extensions#tabline#enabled = 1

Plug 'sheerun/vim-polyglot'

Plug 'dikiaap/minimalist'
Plug 'tomasr/molokai'
Plug 'dracula/vim', { 'as': 'dracula' }

call plug#end()

autocmd VimEnter * if len(filter(values(g:plugs), '!isdirectory(v:val.dir)'))
  \| PlugInstall --sync | source $MYVIMRC
\| endif

cd ~/Desktop/code

syntax on
color dracula

set termguicolors
set nu
set rnu
set expandtab
set tabstop=2
set shiftwidth=2
set mouse=a

inoremap {<CR> {<CR>}<Esc>O

let g:template =<< END
#include <bits/stdc++.h>
using namespace std;

#ifdef LOCAL
#include "templates/debug.h"
#else
#define dbg(x...)
#endif

int main() {
  ios::sync_with_stdio(0), cin.tie(0);
  int tt; cin >> tt;
  for (int tc = 1; tc <= tt; tc++) {
    cout << "Case #" << tc << ": ";
  }
}
END

autocmd filetype cpp nnoremap <C-N> :<C-U>%d \| call setline(1, g:template)<CR>G2k$
autocmd filetype cpp nnoremap <C-B> :<C-U>w \| !g++ -O2 -DLOCAL -std=c++17 %:r.cpp -o run<CR>
autocmd filetype cpp nnoremap <C-C> :<C-U>call writefile(split(getreg('+'), '\n'), 'in') \| !start cmd /c "run < in & pause"<CR><CR>
autocmd filetype cpp nnoremap <C-S> :<C-U>!start cmd /c "run & timeout /t -1 /nobreak"<CR><CR>
nnoremap <C-A> :<C-U>%y+<CR>
