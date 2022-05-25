import React, { useState, useEffect } from 'react';
import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';
import PlusIcon from 'components/general/PlusIcon';
import InstitutionCreationDialog from 'components/institutions/InstitutionCreationDialog';
import InstitutionList from 'components/institutions/InstitutionList';
import InstitutionService from 'shared/services/Institution.service';

function InstitutionOverview() {
  const [institutions, setInstitutions] = useState([]);
  const [isLoading, setIsLoading] = useState(false);
  useEffect(() => {
    setIsLoading(true);
    InstitutionService.getMyInstitutions()
      .then((newInstitutions) => {
        setInstitutions(newInstitutions);
        setIsLoading(false);
      })
      .catch((e) => {
        console.error(e);
        if (e.response?.data) {
          alert(e.response.data);
        }
      })
      .finally(() => setIsLoading(false));
  }, []);

  function onLeft(institution) {
    setInstitutions(institutions.filter((i) => i.id !== institution.id));
  }

  const [createOpen, setCreateOpen] = useState(false);

  return (
    <>
      <Stack direction="row" className="stack" alignItems="Center">
        <Typography variant="h5" sx={{ pr: 1 }}>Your Institutions</Typography>
        <PlusIcon onClick={() => setCreateOpen(true)} />
      </Stack>
      <br />
      <InstitutionList
        isLoading={isLoading}
        institutions={institutions}
        onLeft={(inst) => onLeft(inst)}
      />
      <InstitutionCreationDialog
        open={createOpen}
        handleClose={() => setCreateOpen(false)}
        onCreated={(newInstitution) => setInstitutions([...institutions, newInstitution])}
      />
    </>
  );
}

export default InstitutionOverview;
