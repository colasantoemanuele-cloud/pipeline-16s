# Le risorse usate dal processo, per il log.
#
# Interfaccia:
#   memoria_di_picco()  -> la memoria residente di picco del processo, come
#                          testo ("1234567 kB"), o "non disponibile"
#
# Il valore viene da /proc/self/status (VmHWM), dove il sistema lo espone. Va
# nel log e non negli artefatti: dipende dalla macchina, e un artefatto deve
# restare identico a parita' di ingresso.

memoria_di_picco <- function() {
  stato <- "/proc/self/status"
  if (!file.exists(stato)) return("non disponibile")
  riga <- grep("^VmHWM:", readLines(stato), value = TRUE)
  if (length(riga) == 0L) "non disponibile" else trimws(sub("^VmHWM:", "", riga))
}
