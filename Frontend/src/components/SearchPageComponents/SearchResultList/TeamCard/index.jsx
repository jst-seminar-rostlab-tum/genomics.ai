import React from 'react';

import { Chip } from '@mui/material';

// import Avatars from 'components/Avatars';
import SearchCard from '../SearchCard';
import LabeledLink from '../LabeledLink';

// import CustomButton from 'components/CustomButton';

// Card to display search result for a single team
// eslint-disable-next-line arrow-body-style
const TeamCard = ({ item: team }) => {
  return (
    <SearchCard
      // action={<Button variant="contained">Join</Button>}
      // action={<CustomButton type="primary">Join</CustomButton>}
      title={team.title}
      link={`/sequencer/teams/${team.id}`}
      primary={
        //  <Tag content={team.visibility} variant="primary-default" />
        <Chip label={team.visibility.toLowerCase()} color="primary" size="small" />
      }
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
