# goldfish-fire-repro

密閉空間火災中水缸金魚存活死因機制的 0D 集總參數模型分析。

從一道網路流傳的腦筋急轉彎題出發(「密閉無窗房間猛烈火災把家具燒成灰燼,水缸裡的金魚死活?」),
嘗試轉譯成一個參數化的工程問題,用簡化物理模型探索:熱死 vs 缺氧死,哪個機制主導?在什麼條件下?

**這不是一份已驗證的工程結論,是一份誠實記錄疊代除錯過程的工作草稿。**
詳見 [`article/article.md`](article/article.md) 第4節「限制與未完成事項」。

## 目錄結構

```
.
├── article/
│   └── article.md          # ReScience C 風格總結文件(含摘要/方法/結果/限制)
├── code/
│   ├── fish_ode_model5.py  # 最終版0D模型(牆體材質可調)
│   ├── sweep.py             # 房間體積 x 水缸距離 相圖掃描
│   ├── sensitivity.py       # χ_r / h_conv / kLa 敏感度測試
│   ├── wall_sensitivity.py  # 牆體材質敏感度測試
│   └── timehistory*.py      # 各修正階段的時間歷程分析
├── figures/                 # 所有輸出圖檔(PNG)
└── README.md
```

## 執行環境

```
pip install numpy scipy matplotlib
```

## 已知限制(摘要,完整版見 article.md)

- 僅單一 Python/scipy 實作,**未與任何獨立工具(CFAST/FDS)交叉驗證**
- 未與實測/實驗數據比對
- 多個關鍵係數(h_conv, χ_r, kLa, 魚類代謝參數)為概估值,未經校準
- 0D 單一均溫假設已知會扭曲對流傳熱的空間(距離)依賴性

## 授權

程式碼與文件皆為分析過程記錄,依個人需求自由使用。
