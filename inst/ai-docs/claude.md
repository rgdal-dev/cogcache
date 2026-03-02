---
title: "Claude restart notes"
format: html
---

- Main design framework and philosophy is in inst/design/docs/gdal-warper-audit.md

Guiding tasks: 

1) Integrate into an explicit caching and config mech that mirrors GDAL cache + config but uses duckdb/parquet  
2) the goal is to *reimplement* the full GDAL warp api as something the community can adopt as  much lither body of work that reflects its proper origins 
3) expose it in an R package *first*, because that's where i'm comfy
