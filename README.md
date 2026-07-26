# goldfish-fire-repro

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21610403.svg)](https://doi.org/10.5281/zenodo.21610403)

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
**結論可信度表**、和**驗證紀錄(Validation Log)**。

這個repo真正想回答的問題,其實不是「金魚有沒有死」,而是:

> **我們怎麼知道自己的答案值得相信?**
> **How do we know our answer is trustworthy?**

## 現在的狀態:Level 3(無因次分析)進行到一半,刻意不急著擴充模型

最新進度:把原本二維的相圖(房間體積 × 水缸-火源距離)嘗試用單一無因次群 $\Pi_{heat}$ collapse——
第一次嘗試測出92.4%準確率,但發現是循環論證(用了模擬輸出當分子);推翻重做後,改用純輸入參數,
收斂到 **95.1%**,但公式裡留了一個**刻意不校準**的自由參數(0.3係數,物理意義未知)。
完整過程見 article.md 的 **VL-05**。

下一步優先順序:
1. 把水缸體積、燃料量、牆體材質也納入同一組collapse測試,看能不能解釋剩下的5%、或找出0.3係數對應的機制
2. 找出這個 $\Pi$ 群會在什麼參數範圍失效
3. 只有在0D+無因次分析的機制層級價值被榨乾之後,才進 Level 4(two-zone / CFAST交叉驗證)
4. CFD 和文獻比對排在更後面

**不寫 fish_ode_model6.py。** 這個專案累積洞見的速度,目前比堆疊模型複雜度更有價值。

## 目錄結構

```
.
├── article/
│   └── article.md          # 完整案例研究:思考演化圖/驗證路線圖/可信度表/Validation Log
├── code/
│   ├── fish_ode_model5.py  # 目前最終版0D模型(牆體材質可調,已知有未修正的結構性限制)
│   ├── sweep.py             # 房間體積 x 水缸距離 相圖掃描
│   ├── sensitivity.py       # χ_r / h_conv / kLa 敏感度測試
│   ├── wall_sensitivity.py  # 牆體材質敏感度測試
│   ├── timehistory*.py      # 各修正階段的時間歷程分析(含已知的bug版本,刻意保留)
│   ├── dimensionless.py     # 無因次分析v1(循環論證版,刻意保留)
│   ├── dimensionless_v2.py  # 無因次分析v2(input-only修正版,95.1%)
│   └── collapse_analysis.py # 2D相圖→1D Π群 collapse可視化
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

MIT License,詳見 [`LICENSE`](LICENSE)。
