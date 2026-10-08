# La curva KatharoSeq e la conformita' dei controlli positivi (S11).
#
# Interfaccia:
#   z_robusto(v)                       lo z modificato di ogni valore
#   conformita(livello, fedelta, profondita, min_controlli)
#                                      per controllo: z di fedelta' e profondita',
#                                      esito (conforme, non conforme, non valutabile)
#   adatta_curva(profondita, fedelta, sensibilita)
#                                      la sigmoide allosterica, la bonta', l'AIC
#                                      e la soglia
#   soglia_determinata(profondita, curva)
#                                      se la soglia e' sostenuta da osservazioni
#
# IL MODELLO. La fedelta' f di un controllo positivo e' la frazione delle sue
# letture assegnate al taxon atteso. In funzione di x = log10(profondita'):
#
#     f(x) = x^h / (k' + x^h)
#
# la sigmoide allosterica di KatharoSeq (Minich et al. 2018), con fedelta'
# massima 1. Si adatta nella forma equivalente f(x) = 1 / (1 + (x50 / x)^h),
# con k' = x50^h: x50 e' il log10 della profondita' a cui la fedelta' vale 1/2.
# k' cresce come una potenza (ordini di grandezza diversi da una curva
# all'altra) mentre x50 resta dell'ordine del log10 delle profondita': con k'
# la stima non converge, con x50 si'. La profondita' minima per una fedelta' s e'
#
#     n* = 10^(x50 * (s / (1 - s))^(1 / h))
#
# LA BONTA'. R^2 = 1 - SS_res / SS_tot, con i residui sulla scala della
# fedelta', non pesati, e SS_tot attorno alla fedelta' media dei punti. Per un
# modello non lineare non e' la frazione di varianza spiegata della
# regressione lineare: non c'e' scomposizione SS_tot = SS_mod + SS_res, il
# valore puo' essere negativo (un modello peggiore della media) e non si
# confronta con un test F. Qui dice soltanto quanto la curva riduce l'errore
# quadratico rispetto alla fedelta' media. Accanto si riporta l'errore
# quadratico medio (RMSE), nella stessa unita' della fedelta'.
#
# L'AIC. Il criterio di Akaike della curva (stats::AIC sul modello nls:
# verosimiglianza gaussiana, due parametri della curva piu' la varianza dei
# residui) serve a confrontare la curva aggregata con l'insieme delle curve per
# piastra, adattate sugli stessi punti: l'AIC del modello per piastra e' la
# somma degli AIC delle sue curve, perche' le piastre sono stimate in modo
# indipendente e le verosimiglianze si moltiplicano. Il confronto vale solo a
# parita' di punti: R/11_controls.R li fissa prima di adattare.
#
# L'ADATTAMENTO. nls con l'algoritmo "port" e h, x50 maggiori di zero. I
# valori iniziali vengono dai dati in modo deterministico, dalla forma
# linearizzata logit(f) = h log(x) - h log(x50) con le fedelta' ristrette a
# [0,01; 0,99]; se la retta non da' una pendenza positiva, h = 10 e x50 = la
# mediana di x. Nessun numero casuale: due esecuzioni danno la stessa curva.
#
# LA CONFORMITA'. Un controllo e' confrontato con gli altri controlli dello
# stesso livello di concentrazione (le cellule seminate), nelle altre piastre:
# e' non conforme se la sua fedelta' e' anomalmente bassa per il suo livello,
# o la sua profondita' anomalamente alta o bassa. Anomalo vuol dire uno z
# modificato oltre Z_MAX (Iglewicz e Hoaglin, 1993): z = 0,6745 (v - mediana) /
# MAD, con la deviazione assoluta mediana; con MAD nullo, 1,2533 volte la
# deviazione assoluta media. La profondita' si confronta in log10(n + 1), e
# anche un controllo senza letture, che non ha fedelta', si confronta con il
# suo livello. Un livello con meno di min_controlli controlli non permette il
# confronto: i suoi controlli sono non valutabili, come quelli senza livello. A bassa concentrazione la fedelta' attesa e'
# bassa, e un controllo dominato dai contaminanti come gli altri del suo
# livello e' conforme.

Z_MAX <- 3.5

z_robusto <- function(v) {
  z <- rep(NA_real_, length(v))
  ok <- !is.na(v)
  if (!any(ok)) return(z)
  m <- stats::median(v[ok])
  mad <- stats::median(abs(v[ok] - m))
  if (mad > 0) {
    z[ok] <- 0.6745 * (v[ok] - m) / mad
  } else {
    media <- mean(abs(v[ok] - m))
    z[ok] <- if (media > 0) (v[ok] - m) / (1.253314 * media) else 0
  }
  z
}

conformita <- function(livello, fedelta, profondita, min_controlli) {
  n <- length(livello)
  z_f <- rep(NA_real_, n)
  z_n <- rep(NA_real_, n)
  esito <- rep("non valutabile", n)
  motivo <- rep("senza livello di concentrazione", n)
  for (l in sort(unique(livello[!is.na(livello)]))) {
    gruppo <- which(!is.na(livello) & livello == l)
    if (length(gruppo) < min_controlli) {
      motivo[gruppo] <- sprintf("livello con %d controlli, meno di %d",
                                length(gruppo), min_controlli)
      next
    }
    # Un controllo senza letture non ha fedelta', ma la sua profondita' si
    # confronta comunque con quella del suo livello.
    z_f[gruppo] <- z_robusto(fedelta[gruppo])
    z_n[gruppo] <- z_robusto(log10(profondita[gruppo] + 1))
    for (i in gruppo) {
      cause <- c(
        if (!is.na(z_f[i]) && z_f[i] < -Z_MAX) "fedelta' bassa per il livello",
        if (z_n[i] > Z_MAX) "profondita' alta per il livello",
        if (z_n[i] < -Z_MAX) "profondita' bassa per il livello"
      )
      esito[i] <- if (length(cause)) "non conforme" else "conforme"
      motivo[i] <- paste(cause, collapse = "; ")
    }
  }
  data.frame(z_fedelta = z_f, z_profondita = z_n, esito = esito, motivo = motivo,
             stringsAsFactors = FALSE)
}

adatta_curva <- function(profondita, fedelta, sensibilita) {
  esito <- list(converge = FALSE, punti = length(profondita), h = NA_real_,
                x50 = NA_real_, k = NA_real_, r2 = NA_real_, rmse = NA_real_,
                aic = NA_real_, soglia = NA_real_, motivo = "")
  if (length(profondita) < 3L) {
    esito$motivo <- sprintf("%d punti, ne servono almeno 3 per due parametri",
                            length(profondita))
    return(esito)
  }
  dati <- data.frame(x = log10(profondita), f = fedelta)
  if (stats::var(dati$f) == 0) {
    esito$motivo <- "fedelta' identica in tutti i punti: nessuna curva"
    return(esito)
  }
  fc <- pmin(pmax(dati$f, 0.01), 0.99)
  retta <- stats::lm(stats::qlogis(fc) ~ log(dati$x))
  pendenza <- unname(stats::coef(retta)[2])
  if (is.finite(pendenza) && pendenza > 0) {
    inizio <- list(h = pendenza,
                   x50 = exp(-unname(stats::coef(retta)[1]) / pendenza))
  } else {
    inizio <- list(h = 10, x50 = stats::median(dati$x))
  }
  modello <- tryCatch(
    stats::nls(f ~ 1 / (1 + (x50 / x)^h), data = dati, start = inizio,
               algorithm = "port", lower = c(h = 1e-6, x50 = 1e-6),
               control = stats::nls.control(maxiter = 500)),
    error = function(e) e
  )
  if (inherits(modello, "error")) {
    esito$motivo <- paste("adattamento non convergente:", conditionMessage(modello))
    return(esito)
  }
  stime <- stats::coef(modello)
  residui <- dati$f - stats::fitted(modello)
  esito$converge <- TRUE
  esito$h <- unname(stime["h"])
  esito$x50 <- unname(stime["x50"])
  esito$k <- esito$x50^esito$h
  esito$r2 <- 1 - sum(residui^2) / sum((dati$f - mean(dati$f))^2)
  esito$rmse <- sqrt(mean(residui^2))
  esito$aic <- as.numeric(stats::AIC(modello))
  esito$soglia <- 10^(esito$x50 * (sensibilita / (1 - sensibilita))^(1 / esito$h))
  esito
}

# LA SOGLIA DETERMINATA DAI DATI. Una curva puo' adattarsi bene e lasciare la
# soglia senza sostegno: se fra la meta' della transizione e la soglia non c'e'
# alcuna osservazione, la posizione della soglia viene solo dalla forma del
# modello, e qualunque posizione nel tratto vuoto darebbe un adattamento
# equivalente (succede quando la curva e' un gradino, con h molto grande, in
# un tratto senza controlli). La soglia e'
# determinata se almeno un'osservazione cade fra il punto medio della curva,
# 10^x50, dove la fedelta' vale 1/2, e la soglia stessa, estremi inclusi: la
# soglia e' allora un'interpolazione fra osservazioni dentro la transizione.
# L'altro lato e' garantito dal criterio di validita' che vuole la soglia
# dentro le profondita' osservate.
#
# Un intervallo di confidenza a profilo della soglia non basta: con un
# gradino quasi perfetto la somma dei quadrati minima e' prossima a zero, la
# forma rigida della sigmoide colloca il gradino in base ai valori di fedelta'
# alle due estremita', e l'intervallo risulta stretto proprio dove i dati non
# dicono nulla.

soglia_determinata <- function(profondita, curva) {
  esito <- list(punto_medio = NA_real_, osservazioni = 0L, determinata = FALSE, motivo = "")
  if (!isTRUE(curva$converge)) {
    esito$motivo <- "curva non stimata"
    return(esito)
  }
  esito$punto_medio <- 10^curva$x50
  estremi <- sort(c(esito$punto_medio, curva$soglia))
  esito$osservazioni <- sum(profondita >= estremi[1] & profondita <= estremi[2])
  if (esito$osservazioni > 0L) {
    esito$determinata <- TRUE
  } else {
    sotto <- profondita[profondita < estremi[1]]
    esito$motivo <- sprintf(paste0(
      "soglia non determinata dai dati: nessuna osservazione fra il punto medio della ",
      "curva (%.0f letture) e la soglia (%.0f)%s"),
      esito$punto_medio, curva$soglia,
      if (length(sotto)) sprintf("; la piu' vicina sotto ha %.0f letture", max(sotto)) else "")
  }
  esito
}
