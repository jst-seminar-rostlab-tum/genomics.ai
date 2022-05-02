import React from 'react';

import { Chip } from '@mui/material';
import SearchCard from '../SearchCard';
import Avatars from 'components/Avatars';

// Card to display search result for a single institution
function InstitutionCard({ item: institution }) {
  return (
    <SearchCard
      title={institution.name}
      avatar={institution.logo}
      link={`/sequencer/institutions/${institution.id}`}
      secondary={`updated on ${institution.updated}`}
      tertiary={(
        <>
          <Avatars
            items={institution.members.map(({ name, image }) => ({ src: image, alt: name }))}
          />
          <Chip
            label={`${institution.membersCount} members`}
            variant="outlined"
            size="small"
            sx={{ color: 'text.secondary' }}
          />
          <Chip
            label={`${institution.teamsCount} teams`}
            variant="outlined"
            size="small"
            sx={{ color: 'text.secondary' }}
          />
        </>
      )}
    />
  );
}

export default InstitutionCard;