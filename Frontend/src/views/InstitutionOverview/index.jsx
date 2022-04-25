import React, { useState, useEffect } from 'react';
import {
  Switch, Route, useRouteMatch,
} from 'react-router-dom';
import Button from '@mui/material/Button';
import HeaderView from 'components/general/HeaderView';
import InstitutionPage from 'views/InstitutionPage';
import styles from './institutionOverview.module.css';
import queryMyInstitutions from 'shared/services/mock/institutions';
import InstitutionCreationDialog from 'components/institutions/InstitutionCreationDialog';
import InstitutionList from 'components/institutions/InstitutionList';

function InstitutionOverview({ sidebarShown }) {
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

  const { path } = useRouteMatch();

  return (
    <Switch>
      <Route exact path={`${path}/`}>
    <>
      <HeaderView
        sidebarShown={sidebarShown}
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
    </Route>

    <Route path={`${path}/:id`}>
      <InstitutionPage sidebarShown={sidebarShown} />
    </Route>
    </Switch>
  );
}

export default InstitutionOverview;
