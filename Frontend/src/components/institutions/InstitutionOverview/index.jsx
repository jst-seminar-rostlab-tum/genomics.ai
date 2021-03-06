import React, { useState } from 'react';
import Stack from '@mui/material/Stack';
import Typography from '@mui/material/Typography';
import PlusIcon from 'components/general/PlusIcon';
import InstitutionCreationDialog from 'components/institutions/InstitutionCreationDialog';
import InstitutionList from 'components/institutions/InstitutionList';
import { useInstitutions } from 'shared/context/institutionContext';

function InstitutionOverview() {
  const { institutions, isLoading, setInstitutions } = useInstitutions();

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
