import React from 'react';

import { Chip, Tooltip } from '@mui/material';

// import Avatars from 'components/Avatars';
import SearchCard from '../SearchCard';
import LabeledLink from '../LabeledLink';

import TeamJoinButton from 'components/teams/detail/TeamJoinButton';
import TeamService from 'shared/services/Team.service';

// Card to display search result for a single team
// eslint-disable-next-line arrow-body-style
const TeamCard = ({ item: team, user, onAction }) => {
  const joinTeam = async (teamToJoin) => {
    await TeamService.joinTeam(teamToJoin.id);
    onAction();
  };

  const visibility = team.visibility.toLowerCase();
  let canJoin = false; // target to remove if backend implements rules checks
  let visibilityTooltip;
  const isMember = team.memberIds.includes(user._id); // check if user already joined

  switch (visibility) {
    case 'public':
      visibilityTooltip = 'Anybody can join the project!';
      canJoin = !isMember;
      break;
    case 'private':
      visibilityTooltip = 'Only invited members can join the project!';
      canJoin = team.invitedMemberIds.includes(user._id); // check if user is invited
      break;
    case 'by institution':
      visibilityTooltip = 'Only institution members can join the project!';
      // check if user is not member but part of institution
      canJoin = !isMember
      && team.institution.memberIds.includes(user._id);
      break;
    default:
      visibilityTooltip = 'unknown';
  }

  return (
    <SearchCard
      action={canJoin && <TeamJoinButton team={team} onJoin={joinTeam} />}
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
              content={team.institution.name}
              to={`/sequencer/institutions/${team.institutionId}`}
            />
          )}
        </>
      )}
    />
  );
};

export default TeamCard;
