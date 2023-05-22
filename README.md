# Fast Multipole Method
> Final project for computer astrophysics 111-2

## N-body gravity simulation -- direct computation
> final_project_direct.cpp
2023/5/22
1. 要先裝SDL2才可以視覺化
	- 我用home brew 裝 : brew install sdl2 
  - 然後用這個compile：g++ final_project_direct.cpp -o final_project_direct.out -L/opt/homebrew/Cellar/sdl2/2.26.5/lib -lSDL2 -I/opt/homebrew/Cellar/sdl2/2.26.5/include (-L跟-I後面都放sdl2裝的位置)
	- 跑程式輸入：![image](https://github.com/ej-1211/CA_final_project/assets/91402798/6cc4ac93-8243-4681-a622-c72956bb92e8)
	- 就可以看到：![image](https://github.com/ej-1211/CA_final_project/assets/91402798/42ccab6e-773f-4f8c-8828-8c05df24be2f) 他就會自己跑
2. 功能：
	- 有加softening distance -> 怕他太接近，跟距離平方反比會出現過大的奇怪力
	- 可以選要不要視覺化，或是直接看結果（太多粒子sdl好像跑不動，所以只好直接印出來），之後測試應該是要用很大的N
	- 有算到底跑了多少實際時間
	- 有印出每一dt的error（但目前還有Bug）
4. 問題：
	- 計算Energy有Bug，應該要是Conserve才對，還沒找到
	- 還沒把每個數字的單位搞清楚
	- Periodic Boundary Condition 有點不懂，是從左邊邊界出去就會從右邊邊界進來嗎？這樣能量是不是就會出一點問題
