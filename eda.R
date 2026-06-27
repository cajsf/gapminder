# eda.R
# Gapminder 데이터 탐색적 분석(EDA) - 심화판
# 사용법: Rscript eda.R
# 입력 : data/gapminder.csv
# 출력 : 콘솔 분석 결과 + document/figures/*.png
#
# 설계 원칙:
#  - 단순평균이 아닌 '인구 가중' 통계를 핵심 지표로 사용 (평균적 인간의 경험)
#  - 분포 형태(왜도/첨도)와 그 시간적 변화를 정량화
#  - 국가 간 불평등의 추이(수렴/발산)를 명시적으로 측정
#  - 소득-수명 회귀의 잔차로 '성과 초과/미달' 국가 식별
#  - 역사적 충격(전쟁/기근/질병)을 데이터에서 자동 탐지

# ---- 0. 설정 & 유틸 ----
input_path <- "data/gapminder.csv"
fig_dir    <- "document/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

section <- function(title) {
  cat("\n", strrep("=", 64), "\n", title, "\n", strrep("=", 64), "\n", sep = "")
}
save_png <- function(name, expr, width = 950, height = 620) {
  path <- file.path(fig_dir, name)
  png(path, width = width, height = height, res = 110)
  on.exit(dev.off())
  force(expr)
  cat(sprintf("  [그래프] %s\n", path))
}
# 표본 왜도/첨도 (초과첨도, excess kurtosis)
skewness <- function(x) { x <- x[is.finite(x)]; n <- length(x); m <- mean(x)
  (sum((x - m)^3) / n) / (sum((x - m)^2) / n)^1.5 }
kurtosis_excess <- function(x) { x <- x[is.finite(x)]; n <- length(x); m <- mean(x)
  (sum((x - m)^4) / n) / (sum((x - m)^2) / n)^2 - 3 }
# 비가중 지니계수
gini <- function(x) { x <- sort(x[is.finite(x) & x >= 0]); n <- length(x)
  if (n < 2) return(NA_real_); (2 * sum(seq_len(n) * x)) / (n * sum(x)) - (n + 1) / n }

df <- read.csv(input_path, stringsAsFactors = FALSE, encoding = "UTF-8")
years <- sort(unique(df$year)); y0 <- min(years); y1 <- max(years); span <- y1 - y0
cont_levels <- sort(unique(df$continent))
cont_col <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
                     cont_levels)

# =====================================================================
section("1. 데이터 개요 & 분포 형태 정량화")
# =====================================================================
cat(sprintf("관측치 %d | 국가 %d | 대륙 %d | 연도 %d (%d~%d, %d년 간격)\n",
            nrow(df), length(unique(df$country)), length(cont_levels),
            length(years), y0, y1, years[2] - years[1]))

cat("\n분포 형태 (전체 풀링):\n")
shape <- data.frame(
  variable = c("lifeExp", "gdpPercap", "log10(gdpPercap)", "pop", "log10(pop)"),
  skewness = round(c(skewness(df$lifeExp), skewness(df$gdpPercap),
                     skewness(log10(df$gdpPercap)), skewness(df$pop),
                     skewness(log10(df$pop))), 3),
  exc_kurt = round(c(kurtosis_excess(df$lifeExp), kurtosis_excess(df$gdpPercap),
                     kurtosis_excess(log10(df$gdpPercap)), kurtosis_excess(df$pop),
                     kurtosis_excess(log10(df$pop))), 3))
print(shape, row.names = FALSE)
cat("\n해석: gdpPercap는 강한 우편향(skew>>0)이며 log 변환 시 거의 대칭화됨\n")
cat("      → 이후 분석에서 GDP는 로그 척도 사용이 타당함\n")

# =====================================================================
section("2. 단순평균 vs 인구가중평균 (왜 구분이 중요한가)")
# =====================================================================
# 단순평균은 모든 국가를 동등 취급 → 소국이 과대대표됨.
# 인구가중평균은 '평균적인 한 사람'이 경험하는 값.
cmp <- do.call(rbind, lapply(years, function(yr) {
  d <- df[df$year == yr, ]
  data.frame(year = yr,
    life_simple = round(mean(d$lifeExp), 1),
    life_wtd    = round(weighted.mean(d$lifeExp, d$pop), 1),
    gdp_simple  = round(mean(d$gdpPercap), 0),
    gdp_wtd     = round(weighted.mean(d$gdpPercap, d$pop), 0))
}))
print(cmp, row.names = FALSE)
d0 <- cmp[cmp$year == y0, ]; dN <- cmp[cmp$year == y1, ]
cat(sprintf("\n%d년 기대수명: 단순평균 %.1f vs 인구가중 %.1f (차이 %.1f세)\n",
            y1, dN$life_simple, dN$life_wtd, dN$life_simple - dN$life_wtd))
cat(sprintf("%d년 1인당GDP: 단순평균 %d vs 인구가중 %d\n",
            y1, dN$gdp_simple, dN$gdp_wtd))
cat("해석: 인구가중 기대수명이 단순평균보다 낮음 → 인구 많은 나라(중국·인도 등)가\n")
cat("      상대적으로 낮은 수명대에 위치. 단순평균은 '세계'를 낙관적으로 왜곡함.\n")

# =====================================================================
section("3. 국가 간 불평등의 추이 (수렴 vs 발산)")
# =====================================================================
ineq <- do.call(rbind, lapply(years, function(yr) {
  d <- df[df$year == yr, ]
  data.frame(year = yr,
    life_sd      = round(sd(d$lifeExp), 2),
    life_p90_p10 = round(quantile(d$lifeExp, .9) - quantile(d$lifeExp, .1), 1),
    loggdp_sd    = round(sd(log10(d$gdpPercap)), 3),   # σ-수렴 지표
    gdp_gini     = round(gini(d$gdpPercap), 3),
    gdp_p90_p10  = round(quantile(d$gdpPercap, .9) / quantile(d$gdpPercap, .1), 1))
}))
print(ineq, row.names = FALSE)
cat(sprintf("\n기대수명 표준편차: %.2f(%d) -> %.2f(%d)  [%s]\n",
            ineq$life_sd[1], y0, ineq$life_sd[nrow(ineq)], y1,
            ifelse(ineq$life_sd[nrow(ineq)] < ineq$life_sd[1], "수렴(격차 축소)", "발산")))
cat(sprintf("log GDP 표준편차 : %.3f(%d) -> %.3f(%d)  [%s]\n",
            ineq$loggdp_sd[1], y0, ineq$loggdp_sd[nrow(ineq)], y1,
            ifelse(ineq$loggdp_sd[nrow(ineq)] < ineq$loggdp_sd[1], "수렴", "발산(격차 확대)")))
cat("해석: 기대수명은 수렴하나 소득(log) 격차는 줄지 않음 → 건강은 따라잡아도\n")
cat("      소득 격차는 고착. '건강 수렴, 소득 발산'의 비대칭.\n")

# =====================================================================
section("4. β-수렴 분석 (가난한 나라가 빨리 성장했는가?)")
# =====================================================================
# 국가별 1952->2007 GDP 연평균성장률(CAGR)을 초기 소득에 회귀.
# 계수<0 이면 초기 빈국이 더 빨리 성장 = 절대수렴.
g <- merge(
  df[df$year == y0, c("country", "continent", "gdpPercap")],
  df[df$year == y1, c("country", "gdpPercap")],
  by = "country", suffixes = c("_0", "_1"))
g$cagr <- (g$gdpPercap_1 / g$gdpPercap_0)^(1 / span) - 1
g$log_gdp0 <- log10(g$gdpPercap_0)
bfit <- lm(cagr ~ log_gdp0, data = g)
beta <- coef(bfit)["log_gdp0"]
cat(sprintf("회귀: CAGR ~ log10(초기GDP)\n  기울기 = %.4f (p = %.3f), R^2 = %.3f\n",
            beta, summary(bfit)$coefficients["log_gdp0", 4], summary(bfit)$r.squared))
cat(sprintf("  결론: 기울기 %s 0 → %s\n",
            ifelse(beta < 0, "<", ">="),
            ifelse(beta < 0, "약한 절대수렴 신호(빈국이 다소 빨리 성장)",
                   "절대수렴 없음(부국이 더 빨리 성장, 발산)")))
cat("\n성장 챔피언 (CAGR 상위 5):\n")
print(head(g[order(-g$cagr), c("country", "continent", "gdpPercap_0", "gdpPercap_1", "cagr")], 5), row.names = FALSE)
cat("\n성장 정체/후퇴 (CAGR 하위 5):\n")
print(head(g[order(g$cagr), c("country", "continent", "gdpPercap_0", "gdpPercap_1", "cagr")], 5), row.names = FALSE)

# =====================================================================
section("5. 소득-수명 관계의 시간적 변화 & 회귀 잔차")
# =====================================================================
cat("연도별 상관 corr(lifeExp, log10 gdp):\n")
cor_by_yr <- sapply(years, function(yr) {
  d <- df[df$year == yr, ]; cor(d$lifeExp, log10(d$gdpPercap)) })
print(round(setNames(cor_by_yr, years), 3))

# 2007년 회귀: 소득으로 설명되는 수명 + 잔차(=성과 초과/미달)
d07 <- df[df$year == y1, ]
fit07 <- lm(lifeExp ~ log10(gdpPercap), data = d07)
d07$resid <- residuals(fit07)
cat(sprintf("\n%d년 회귀 lifeExp ~ log10(gdp): R^2 = %.3f\n", y1, summary(fit07)$r.squared))
cat("소득 대비 '기대수명 초과 성과' 상위 5 (잔차 +):\n")
print(head(d07[order(-d07$resid), c("country", "continent", "gdpPercap", "lifeExp", "resid")], 5), row.names = FALSE)
cat("소득 대비 '기대수명 미달' 상위 5 (잔차 -):\n")
print(head(d07[order(d07$resid), c("country", "continent", "gdpPercap", "lifeExp", "resid")], 5), row.names = FALSE)
cat("해석: 잔차 +는 '저소득에도 수명이 높은' 보건 효율 국가, -는 그 반대.\n")

# =====================================================================
section("6. 역사적 충격 자동 탐지 (연속 시점 간 급락)")
# =====================================================================
ord <- df[order(df$country, df$year), ]
ord$d_life <- ave(ord$lifeExp,  ord$country, FUN = function(x) c(NA, diff(x)))
ord$d_gdp_pct <- ave(ord$gdpPercap, ord$country,
                     FUN = function(x) c(NA, diff(x) / head(x, -1) * 100))
cat("기대수명 최대 급락 7건 (5년 사이 하락폭):\n")
top_life_drop <- head(ord[order(ord$d_life), c("country", "year", "lifeExp", "d_life")], 7)
print(top_life_drop, row.names = FALSE)
cat("\n1인당 GDP 최대 급락 7건 (%):\n")
top_gdp_drop <- head(ord[order(ord$d_gdp_pct), c("country", "year", "gdpPercap", "d_gdp_pct")], 7)
top_gdp_drop$d_gdp_pct <- round(top_gdp_drop$d_gdp_pct, 1)
print(top_gdp_drop, row.names = FALSE)
cat("해석: 르완다('92→'97 학살), 캄보디아(크메르루주), 남부아프리카(HIV/AIDS),\n")
cat("      쿠웨이트(걸프전 유가/생산 붕괴) 등 실제 역사적 사건과 일치.\n")

# =====================================================================
section("7. 대륙 내 이질성 (평균이 숨기는 분산)")
# =====================================================================
het <- do.call(rbind, lapply(cont_levels, function(cn) {
  d <- d07[d07$continent == cn, ]
  data.frame(continent = cn, n = nrow(d),
    life_min = round(min(d$lifeExp), 1), life_med = round(median(d$lifeExp), 1),
    life_max = round(max(d$lifeExp), 1),
    gdp_ratio_max_min = round(max(d$gdpPercap) / min(d$gdpPercap), 1))
}))
cat(sprintf("%d년 대륙 내 분포:\n", y1))
print(het, row.names = FALSE)
cat("해석: 아시아는 내부 격차가 극심(최부국/최빈국 GDP 비율 최대) → '아시아 평균'은\n")
cat("      일본·쿠웨이트와 아프가니스탄·미얀마를 뭉뚱그린 무의미한 수치.\n")

# =====================================================================
section("8. 그래프 생성")
# =====================================================================

# 8-1. 기대수명 분포의 시간적 변화 (이봉 -> 단봉 수렴) : 연도별 밀도 중첩
save_png("01_lifeExp_density_by_year.png", {
  sel <- c(1952, 1977, 2007)
  cols <- c("#bdbdbd", "#fb8072", "#1f78b4")
  dens <- lapply(sel, function(yr) density(df$lifeExp[df$year == yr]))
  plot(NA, xlim = c(20, 90), ylim = c(0, max(sapply(dens, function(d) max(d$y)))),
       xlab = "Life expectancy (years)", ylab = "Density",
       main = "Life expectancy distribution: bimodal -> unimodal")
  for (i in seq_along(sel)) { polygon(dens[[i]], col = adjustcolor(cols[i], .4), border = cols[i], lwd = 2) }
  legend("topleft", legend = sel, fill = adjustcolor(cols, .5), border = cols, bty = "n")
})

# 8-2. 단순평균 vs 인구가중평균 추세
save_png("02_simple_vs_weighted.png", {
  par(mfrow = c(1, 2))
  matplot(cmp$year, cmp[, c("life_simple", "life_wtd")], type = "b", pch = 16,
          lty = 1, col = c("#999999", "#d62728"), lwd = 2,
          xlab = "Year", ylab = "Life expectancy", main = "Life exp: simple vs pop-weighted")
  legend("bottomright", c("Simple mean", "Pop-weighted"),
         col = c("#999999", "#d62728"), lwd = 2, pch = 16, bty = "n")
  matplot(cmp$year, cmp[, c("gdp_simple", "gdp_wtd")], type = "b", pch = 16,
          lty = 1, col = c("#999999", "#2ca02c"), lwd = 2,
          xlab = "Year", ylab = "GDP per capita", main = "GDP: simple vs pop-weighted")
  legend("topleft", c("Simple mean", "Pop-weighted"),
         col = c("#999999", "#2ca02c"), lwd = 2, pch = 16, bty = "n")
  par(mfrow = c(1, 1))
}, width = 1200)

# 8-3. 불평등 추이 (수명 SD 수렴 vs log GDP SD 정체)
save_png("03_inequality_trend.png", {
  par(mar = c(5, 4, 4, 5))
  plot(ineq$year, ineq$life_sd, type = "b", pch = 16, col = "#d62728", lwd = 2,
       xlab = "Year", ylab = "SD of life expectancy", main = "Inequality: health converges, income does not")
  par(new = TRUE)
  plot(ineq$year, ineq$loggdp_sd, type = "b", pch = 17, col = "#1f78b4", lwd = 2,
       axes = FALSE, xlab = "", ylab = "")
  axis(4); mtext("SD of log10(GDP)", side = 4, line = 3)
  legend("top", c("Life exp SD (left)", "log GDP SD (right)"),
         col = c("#d62728", "#1f78b4"), lwd = 2, pch = c(16, 17), bty = "n")
})

# 8-4. β-수렴 산점도
save_png("04_beta_convergence.png", {
  plot(g$log_gdp0, g$cagr * 100, col = adjustcolor(cont_col[g$continent], .8), pch = 19,
       xlab = sprintf("log10(GDP per capita in %d)", y0),
       ylab = sprintf("Annualized GDP growth %d-%d (%%)", y0, y1),
       main = "Beta-convergence: initial income vs growth")
  abline(lm(I(cagr * 100) ~ log_gdp0, data = g), col = "grey30", lty = 2, lwd = 2)
  legend("topright", legend = cont_levels, col = cont_col[cont_levels], pch = 19, bty = "n")
})

# 8-5. 소득-수명 회귀 잔차 (성과 초과/미달, 2007)
save_png("05_life_residuals_2007.png", {
  o <- d07[order(d07$resid), ]
  ext <- rbind(head(o, 8), tail(o, 8))
  cols <- ifelse(ext$resid >= 0, "#2ca02c", "#d62728")
  par(mar = c(5, 9, 4, 2))
  barplot(ext$resid, horiz = TRUE, names.arg = ext$country, las = 1, col = cols,
          cex.names = 0.8, xlab = "Residual: actual - predicted life exp (years)",
          main = sprintf("Health performance vs income (%d)", y1))
  abline(v = 0, col = "grey40")
})

# 8-6. 상관 변화 + 산점도(2007, 인구 가중 크기)
save_png("06_life_vs_gdp_2007.png", {
  plot(d07$gdpPercap, d07$lifeExp, log = "x",
       col = adjustcolor(cont_col[d07$continent], .7), pch = 19,
       cex = sqrt(d07$pop) / 8000,
       xlab = "GDP per capita (log scale)", ylab = "Life expectancy",
       main = sprintf("Life exp vs GDP (%d); point size = population", y1))
  curve(coef(fit07)[1] + coef(fit07)[2] * log10(x), add = TRUE, col = "grey30", lty = 2, lwd = 2)
  legend("bottomright", legend = cont_levels, col = cont_col[cont_levels], pch = 19, bty = "n")
})

# 8-7. 주요 충격 국가의 기대수명 궤적
save_png("07_shock_trajectories.png", {
  watch <- c("Rwanda", "Cambodia", "China", "Zimbabwe", "Botswana")
  watch <- watch[watch %in% df$country]
  pal <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")
  plot(NA, xlim = range(years), ylim = c(20, 85),
       xlab = "Year", ylab = "Life expectancy",
       main = "Trajectories through historical shocks")
  for (i in seq_along(watch)) {
    s <- df[df$country == watch[i], ]; s <- s[order(s$year), ]
    lines(s$year, s$lifeExp, col = pal[i], lwd = 2, type = "b", pch = 16)
  }
  legend("bottomleft", legend = watch, col = pal[seq_along(watch)], lwd = 2, pch = 16, bty = "n")
})

# 8-8. 대륙 내 이질성 (2007 기대수명 박스 + 개별점)
save_png("08_within_continent_spread.png", {
  boxplot(lifeExp ~ continent, data = d07, col = adjustcolor(cont_col[cont_levels], .4),
          outline = FALSE, ylab = "Life expectancy",
          main = sprintf("Within-continent spread (%d)", y1))
  set.seed <- NULL
  for (i in seq_along(cont_levels)) {
    yv <- d07$lifeExp[d07$continent == cont_levels[i]]
    xv <- rep(i, length(yv)) + seq(-0.15, 0.15, length.out = length(yv))
    points(xv, yv, col = cont_col[cont_levels[i]], pch = 19, cex = 0.7)
  }
})

section("EDA 완료")
cat(sprintf("그래프 8종이 '%s'에 저장되었습니다.\n", fig_dir))
