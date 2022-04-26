import React, { useState, useEffect } from 'react';
import Button from '@mui/material/Button';
import HeaderView from 'components/general/HeaderView';
import InstitutionCard from 'components/institutionOverview/InstitutionCard';
import styles from './institutionOverview.module.css';
import queryMyInstitutions from 'shared/services/mock/institutions';
import InstitutionCreationDialog from 'components/institutionOverview/InstitutionCreationDialog';

function InstitutionOverview({ sidebarShown }) {
  const [institutions, setInstitutions] = useState([]);
  useEffect(() => {
    queryMyInstitutions()
      .then((newInstitutions) => setInstitutions(newInstitutions))
      .catch((ignored) => { console.log(ignored); });
  }, [setInstitutions]);

  function onLeft(institution) {
    setInstitutions(institutions.filter((i) => i.id !== institution.id));
  }

  const [createOpen, setCreateOpen] = useState(false);

  return (
    <>
      <HeaderView
        sidebarShown={sidebarShown}
        title="My Institutions"
        replaceHeaderRight={(
          <Button onClick={() => setCreateOpen(true)}>Create</Button>
        )}
      >
        <div className={styles.content}>
          {institutions.length === 0 ? 'No institutions.' : ''}
          {institutions.map((institution) => (
            <div key={institution.id}>
              <InstitutionCard
                institution={institution}
                onLeft={(inst) => onLeft(inst)}
              />
              <div className={styles.cardSpacing} />
            </div>
          ))}
        </div>
      </HeaderView>
      <InstitutionCreationDialog
        open={createOpen}
        handleClose={() => setCreateOpen(false)}
        onCreated={(newInstitution) => setInstitutions([...institutions, newInstitution])}
      />
    </>
  );
}

export default InstitutionOverview;
