# Configurazione Collaborazione - Istruzioni Finali

Questo documento contiene le istruzioni per completare la configurazione della strategia di collaborazione su GitHub.

## Stato Attuale

Ho completato le seguenti operazioni:

- ✅ Creato `CONTRIBUTING.md` con linee guida dettagliate per i collaboratori
- ✅ Aggiornato `README.md` con la strategia di branching
- ✅ Creato il branch `develop` localmente
- ✅ Fatto commit e push delle modifiche

## Passaggi da Completare su GitHub

### 1. Creare il Branch `develop` su GitHub

Il branch `develop` è stato creato localmente ma deve essere pubblicato su GitHub. Esegui questi comandi:

```bash
# Assicurati di essere sul branch develop
git checkout develop

# Aggiorna con l'ultima versione di main
git pull origin main

# Pusha il branch develop (da fare con le tue credenziali)
git push origin develop
```

**Alternativa via GitHub Web**:
1. Vai su https://github.com/gbucci/koios_vrs
2. Clicca sul menu a tendina dei branch (dove c'è scritto "main")
3. Digita "develop" nel campo di ricerca
4. Clicca su "Create branch: develop from main"

### 2. Proteggere il Branch `main`

Per impedire push diretti al branch `main` e preservare le versioni su Zenodo:

1. Vai su **Settings** del repository
   - URL: https://github.com/gbucci/koios_vrs/settings

2. Nel menu laterale, clicca su **Branches**

3. Nella sezione **Branch protection rules**, clicca su **Add rule**

4. Configura la regola:
   - **Branch name pattern**: `main`
   - Seleziona le seguenti opzioni:
     - ✅ **Require a pull request before merging**
       - ✅ Require approvals (almeno 1)
     - ✅ **Require status checks to pass before merging** (se hai CI/CD)
     - ✅ **Include administrators** (importante!)
     - ✅ **Do not allow bypassing the above settings**

5. Clicca su **Create** per salvare

### 3. Impostare il Branch di Default su `develop`

Per facilitare il lavoro dei collaboratori, imposta `develop` come branch di default:

1. Vai su **Settings** > **Branches**
2. Nella sezione **Default branch**, clicca sull'icona di cambio (frecce)
3. Seleziona `develop` dal menu a tendina
4. Clicca su **Update**
5. Conferma l'operazione

### 4. Aggiungere Collaboratori

Per aggiungere collaboratori al repository:

1. Vai su **Settings** > **Collaborators and teams**
   - URL: https://github.com/gbucci/koios_vrs/settings/access

2. Clicca su **Add people**

3. Inserisci il nome utente GitHub o email del collaboratore

4. Seleziona il livello di accesso:
   - **Write**: Può pushare su branch non protetti, creare PR
   - **Maintain**: Write + gestione issues e settings base
   - **Admin**: Accesso completo (sconsigliato per collaboratori normali)

5. Per la maggior parte dei collaboratori, scegli **Write**

### 5. Creare una Pull Request per Integrare le Modifiche

Le modifiche sono state pushate sul branch `claude/setup-collaboration-access-01JyfMPVM9mAFNLUmsnP3wdS`. Per integrarle:

**Opzione A - Merge in `develop`** (consigliato):
1. Vai su https://github.com/gbucci/koios_vrs/pull/new/claude/setup-collaboration-access-01JyfMPVM9mAFNLUmsnP3wdS
2. Imposta **base: develop** (non main!)
3. Crea la PR e mergela

**Opzione B - Merge in `main`** (se vuoi includerlo nella prossima release):
1. Vai su https://github.com/gbucci/koios_vrs/pull/new/claude/setup-collaboration-access-01JyfMPVM9mAFNLUmsnP3wdS
2. Imposta **base: main**
3. Crea la PR e mergela
4. Successivamente, aggiorna `develop` da `main`:
   ```bash
   git checkout develop
   git merge main
   git push origin develop
   ```

## Workflow di Collaborazione Finale

Una volta completati tutti i passaggi, il workflow sarà:

```
main (protetto, solo release Zenodo)
  ↑
  | PR per release
  |
develop (sviluppo principale)
  ↑
  | PR da feature branches
  |
feature/nuova-funzionalita (lavoro dei collaboratori)
```

## Istruzioni per i Collaboratori

Condividi queste istruzioni con i tuoi collaboratori:

1. **Clone del repository**:
   ```bash
   git clone https://github.com/gbucci/koios_vrs.git
   cd koios_vrs
   ```

2. **Checkout del branch develop**:
   ```bash
   git checkout develop
   git pull origin develop
   ```

3. **Creazione di un feature branch**:
   ```bash
   git checkout -b feature/nome-funzionalita
   ```

4. **Sviluppo e commit**:
   ```bash
   # ... lavoro sul codice ...
   git add .
   git commit -m "Descrizione modifiche"
   ```

5. **Push e Pull Request**:
   ```bash
   git push origin feature/nome-funzionalita
   # Poi creare PR su GitHub verso 'develop'
   ```

## Rilascio di Nuove Versioni su Zenodo

Quando sei pronto per una nuova release:

1. Testa tutto su `develop`
2. Crea una PR da `develop` a `main`
3. Mergia la PR (dopo review se hai abilitato le protezioni)
4. Tagga la release su `main`:
   ```bash
   git checkout main
   git pull origin main
   git tag -a v1.1.0 -m "Descrizione release"
   git push origin v1.1.0
   ```
5. Crea la release su GitHub
6. Zenodo rileverà automaticamente il nuovo tag

## Verifica della Configurazione

Per verificare che tutto sia configurato correttamente:

1. ✅ Il branch `develop` esiste su GitHub
2. ✅ Il branch `main` è protetto (vai su Settings > Branches per verificare)
3. ✅ Il branch di default è `develop`
4. ✅ I collaboratori sono stati aggiunti
5. ✅ Le modifiche (CONTRIBUTING.md, README.md) sono state mergiate

## Supporto

Per qualsiasi problema con la configurazione, consulta:
- [GitHub Docs - Branch Protection](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-protected-branches/about-protected-branches)
- [GitHub Docs - Managing Access](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/managing-repository-settings/managing-teams-and-people-with-access-to-your-repository)
