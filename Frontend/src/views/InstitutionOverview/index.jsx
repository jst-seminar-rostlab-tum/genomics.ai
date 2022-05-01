import React, { useState, useEffect } from 'react';
import Button from '@mui/material/Button';
import HeaderView from 'components/general/HeaderView';
import styles from './institutionOverview.module.css';
import queryMyInstitutions from 'shared/services/mock/institutions';
import InstitutionCreationDialog from 'components/institutions/InstitutionCreationDialog';
import InstitutionList from 'components/institutions/InstitutionList';

function InstitutionOverview() {
  const [institutions, setInstitutions] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  useEffect(async () => {
    setIsLoading(true);
    const newInstitutions = await queryMyInstitutions();
    setInstitutions(newInstitutions);
    setIsLoading(false);
  }, [setInstitutions, setIsLoading]);

  function onLeft(institution) {
    setInstitutions(institutions.filter((i) => i.id !== institution.id));
  }

  const [createOpen, setCreateOpen] = useState(false);

  return (
    <>
      <HeaderView
        title="My Institutions"
        replaceHeaderRight={(
          <Button onClick={() => setCreateOpen(true)}>Create</Button>
        )}
      >
        <div className={styles.content}>
          <InstitutionList
            isLoading={isLoading}
            institutions={institutions}
            onLeft={(inst) => onLeft(inst)}
          />
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
