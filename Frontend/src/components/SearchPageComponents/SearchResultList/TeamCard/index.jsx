import React from 'react';

import { Chip } from '@mui/material';

// import Avatars from 'components/Avatars';
import SearchCard from '../SearchCard';
import LabeledLink from '../LabeledLink';

import TeamJoinButton from 'components/teams/detail/TeamJoinButton';
import TeamService from 'shared/services/Team.service';

import { useAuth } from 'shared/context/authContext';

// Card to display search result for a single team
// eslint-disable-next-line arrow-body-style
const TeamCard = ({ item: team }) => {
  const joinTeam = (team) => {
    TeamService.joinTeam(team.id);
  };

  const [user] = useAuth();
  const isMember = team.memberIds.includes(user._id);

  return (
    <SearchCard
      action={!isMember && <TeamJoinButton team={team} onJoin={joinTeam} />}
      title={team.title}
      link={`/sequencer/teams/${team.id}`}
      primary={
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
