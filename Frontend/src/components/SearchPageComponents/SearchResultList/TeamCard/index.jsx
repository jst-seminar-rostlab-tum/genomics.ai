import React from 'react';

import { Chip, Tooltip } from '@mui/material';

// import Avatars from 'components/Avatars';
import SearchCard from '../SearchCard';
import LabeledLink from '../LabeledLink';

// import CustomButton from 'components/CustomButton';

// Card to display search result for a single team
// eslint-disable-next-line arrow-body-style
const TeamCard = ({ item: team }) => {
  const visibility = team.visibility.toLowerCase();
  let visibilityTooltip;

  switch (visibility) {
    case 'public':
      visibilityTooltip = 'Anybody can join the project!';
      break;
    case 'private':
      visibilityTooltip = 'Only invited members can join the project!';
      break;
    case 'by institution':
      visibilityTooltip = 'Only institution members can join the project!';
      break;
    default:
      visibilityTooltip = 'unknown';
  }

  return (
    <SearchCard
      // action={<CustomButton type="primary">Join</CustomButton>}
      title={team.title}
      link={`/sequencer/teams/${team.id}`}
      primary={(
        // <Tag content={team.visibility} variant="primary-default" />
        <Tooltip title={visibilityTooltip} placement="right-end">
          <Chip label={visibility} color="primary" size="small" />
        </Tooltip>
        )}
      // secondary={`updated on ${team.updated}`}
      tertiary={(
        <>
          {/* <Avatars
            items={team.members.map(({ name, image }) => ({ src: image, alt: name }))}
          /> */}
          <Chip
            label={`${team.memberIds.length + team.adminIds.length} members`}
            variant="outlined"
            size="small"
            sx={{ color: 'text.secondary' }}
          />
          <Chip
            label={`${team.projects.length} projects`}
            variant="outlined"
            size="small"
            sx={{ color: 'text.secondary' }}
          />
          {team.institutionId && (
            <LabeledLink
              label="Institution"
              tooltip="The institution managing the team."
              content={team.institutionTitle}
              to={`/sequencer/institutions/${team.institutionId}`}
            />
          )}
        </>
      )}
    />
  );
};

export default TeamCard;
