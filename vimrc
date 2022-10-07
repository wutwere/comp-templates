call plug#begin()

Plug 'dense-analysis/ale'
let g:ale_linters = {'cpp': ['g++']}
let g:ale_cpp_cc_executable = 'g++'
let g:ale_cpp_cc_options = '-std=c++17 -Wall -Wextra -DLOCAL'

Plug 'vim-airline/vim-airline'
let g:airline_powerline_fonts = 1
let g:airline_theme = 'edge'
let g:airline#extensions#tabline#enabled = 1

Plug 'sheerun/vim-polyglot'

Plug 'sainnhe/edge'
Plug 'arzg/vim-colors-xcode'
Plug 'tomasr/molokai'
Plug 'joshdick/onedark.vim'
Plug 'dracula/vim', { 'as': 'dracula' }

call plug#end()

autocmd VimEnter * if len(filter(values(g:plugs), '!isdirectory(v:val.dir)'))
  \| PlugInstall --sync | source $MYVIMRC
\| endif

cd ~/Desktop/code

syntax on
color molokai

set termguicolors
set nu
set rnu
set expandtab
set tabstop=2
set shiftwidth=2

inoremap {<CR> {<CR>}<Esc>O

autocmd filetype cpp nnoremap <C-N> ggdGa#include <bits/stdc++.h><CR>using namespace std;<CR><CR>#ifdef LOCAL<CR>#include "templates/debug.h"<CR>#else<CR>#define dbg(x...)<CR>#endif<CR><CR>int main() {<CR>ios::sync_with_stdio(0), cin.tie(0);<CR>}<Esc>Oint tt; cin >> tt;<CR>for (int tc = 1; tc <= tt; tc++) {<CR>}<Esc>Ocout << "Case #" << tc << ": ";<Esc>
autocmd filetype cpp nnoremap <C-B> :<C-U>w <bar> !g++ -O2 -DLOCAL -std=c++17 %:r.cpp -o run<CR>
autocmd filetype cpp nnoremap <C-C> :<C-U>1 split in<CR>VG$"+p:w<CR>:bd<CR>:!start cmd /c "run < in & pause"<CR><CR>
autocmd filetype cpp nnoremap <C-S> :<C-U>!start cmd /c "run & pause"<CR><CR>
autocmd filetype cpp nnoremap <C-A> :<C-U>:%y+<CR>
