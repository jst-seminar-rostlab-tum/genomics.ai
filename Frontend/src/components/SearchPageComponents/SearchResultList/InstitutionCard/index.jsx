import React from 'react';

import { Chip } from '@mui/material';
import SearchCard from '../SearchCard';
import { formatDate } from 'shared/utils/common/utils';

// Card to display search result for a single institution
function InstitutionCard({ item: institution }) {
  return (
    <SearchCard
      title={institution.name}
      avatar={institution.profilePictureURL}
      displayAvatar
      link={`/sequencer/institutions/${institution.id}`}
      secondary={institution.updatedAt && `updated on ${formatDate(institution.updatedAt)}`}
      tertiary={(
        <>
          <Chip
            label={`${institution.memberIds.length} members`}
            variant="outlined"
            size="small"
            sx={{ color: 'text.secondary' }}
          />
          <Chip
            label={`${institution.teams.length} teams`}
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
