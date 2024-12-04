call plug#begin()

Plug 'itchyny/lightline.vim'
let g:lightline = { 'colorscheme': 'molokai' }

Plug 'sheerun/vim-polyglot'

Plug 'cormacrelf/vim-colors-github'
Plug 'tomasr/molokai'
Plug 'dracula/vim', { 'as': 'dracula' }

Plug 'ThePrimeagen/vim-be-good'

call plug#end()

autocmd VimEnter * if len(filter(values(g:plugs), '!isdirectory(v:val.dir)'))
  \| PlugInstall --sync | source $MYVIMRC
\| endif

cd ~/dev/comp

syntax on
color molokai

set background=dark
set termguicolors
set nu
set rnu
set expandtab
set tabstop=2
set shiftwidth=2
set mouse=a
set cindent cino=j1,(0,ws,Ws

inoremap {<CR> {<CR>}<Esc>O

let g:template =<< END
#include <bits/stdc++.h>
using namespace std;

#ifdef LOCAL
#include "templates/debug.h"
#else
#define dbg(x...)
#endif

#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define all(x) begin(x), end(x)
#define sz(x) (int)(x).size()
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;

int main() {
  cin.tie(0)->sync_with_stdio(0);
  cin.exceptions(cin.failbit);
  for (int __, _ = (cin >> __, 0); ++_ <= __;) {
    cout << "Case #" << _ << ": ";
  }
}
END

function Run_clipboard()
  call writefile(split(getreg('+'), '\n'), 'in.txt')
  if expand('%:e') == 'cpp'
    execute '!start cmd /c "run < in.txt & pause"'
  elseif expand('%:e') == 'py'
    execute 'w | !start cmd /c "py ' . expand('%:t') . ' < in.txt & pause"'
  endif
endfunction

autocmd filetype cpp nnoremap <C-N> :%d \| call setline(1, g:template)<CR>G2k
autocmd filetype cpp nnoremap <C-B> gg=G'':w \| !g++ -O2 -DLOCAL -std=c++17 %:r.cpp -o run<CR>
nnoremap <C-C> :call Run_clipboard()<CR><CR>
autocmd filetype cpp nnoremap <C-S> :!start cmd /c "run & timeout /t -1 /nobreak"<CR><CR>
nnoremap <C-A> :%y+<CR>
