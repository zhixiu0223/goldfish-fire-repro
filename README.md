# goldfish-fire-repro

**Computational Thinking Case Study — not a fire model of goldfish survival.**

From a popular brain teaser ("Can goldfish survive a severe fire in a sealed room?"),
this repository demonstrates how an ill-posed question can be translated into a
parameterized engineering problem, and how model assumptions, sensitivity analysis,
and validation progressively reshape the conclusions — including the several times
the model turned out to be wrong, and why.

中文版定位:這不是一個「金魚火災存活研究」,而是一份**工程建模與驗證方法論的案例研究**——
記錄一道資訊不足、答案本身無單一正解的腦筋急轉彎,如何被轉譯成可參數化、可反駁、可逐步修正
的工程模型,以及過程中模型被證明「不可信」的每一個時刻。

👉 **如果你只想知道「金魚死不死」,這裡沒有確定答案,而且是刻意沒有。**
誠實的結論是:在目前的驗證程度下,還不到能給確定答案的階段。真正值得看的是
[`article/article.md`](article/article.md) 裡的**思考演化圖**、**驗證路線圖**、
**結論可信度表**、和**失敗紀錄(Failure Log)**。

## 現在的狀態:刻意停在整理階段,不再擴充模型

上一輪的產出不是「更複雜的模型」,而是回頭做這幾件事:
1. 整理哪些結論穩健、哪些不穩健(見 article.md 的可信度表)
2. 把每一次模型被推翻的過程完整記錄下來(見 Failure Log)
3. 畫出目前在 Validation Roadmap 上的實際位置(Level 2 完成,Level 3 未做)

**下一步是無因次分析(Level 3),不是 fish_ode_model6.py。**

## 目錄結構

```
.
├── article/
│   └── article.md          # 完整案例研究:思考演化圖/驗證路線圖/可信度表/失敗紀錄
├── code/
│   ├── fish_ode_model5.py  # 目前最終版0D模型(牆體材質可調,已知有未修正的結構性限制)
│   ├── sweep.py             # 房間體積 x 水缸距離 相圖掃描
│   ├── sensitivity.py       # χ_r / h_conv / kLa 敏感度測試
│   ├── wall_sensitivity.py  # 牆體材質敏感度測試
│   └── timehistory*.py      # 各修正階段的時間歷程分析(含已知的bug版本,刻意保留)
├── figures/                 # 所有輸出圖檔(PNG)
└── README.md
```

## 執行環境

```
pip install numpy scipy matplotlib
```

## 已知限制(完整版見 article.md)

- 僅單一 Python/scipy 實作,**未與任何獨立工具(NIST CFAST/FDS)交叉驗證**
- 未與實測/實驗數據比對
- 多個關鍵係數(h_conv, χ_r, kLa, 魚類代謝參數)為概估值,未經校準
- 0D 單一均溫假設已知會扭曲對流傳熱的空間(距離)依賴性,水缸-火源距離的定量結果不可信
- 半無限固體牆體傳導假設,只在特定時間窗內成立(單層石膏板僅2.4分鐘),超過此窗口目前是權宜處理

## 授權

程式碼與文件皆為分析過程記錄,依個人需求自由使用。
