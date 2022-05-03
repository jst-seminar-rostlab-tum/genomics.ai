import * as React from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import {
  Box, FormControl, InputLabel, MenuItem, Modal, Select, Stack,
} from '@mui/material';
import CircleIcon from '@mui/icons-material/Circle';
import { useHistory } from 'react-router-dom';
import DeleteOutlineIcon from '@mui/icons-material/DeleteOutline';

export default function ProjectBarCard({ projectId, name, status }) {
  const history = useHistory();
  const [color, setColor] = React.useState(status === 'DONE' ? 'lightGreen' : status === 'IN PROGRESS' ? 'orange' : 'red');
  const [typographyColor, setTypographyColor] = React.useState(status === 'UPLOAD FAILED' ? 'red' : 'black');
  const [addTeam, setAddTeam] = React.useState(false);
  const [team, setTeam] = React.useState('');
  const handleOpen = () => setAddTeam(true);
  const handleClose = () => setAddTeam(false);
  const handleChange = (event) => {
    setTeam(event.target.value);
    window.localStorage.setItem('Team', event.target.value);
  };

  React.useEffect(() => {
    if (window.localStorage.getItem('Team')) { setTeam(window.localStorage.getItem('Team')); }
  }, []);

  return (
    <Card sx={{
      marginTop: '5', marginBottom: '0.5em', borderStyle: 'solid', borderColor: '#C8C8C8', borderWidth: '0.1px',
    }}
    >
      <Stack direction="row" sx={{ justifyContent: 'space-between', height: '56px' }}>
        <Stack
          direction="row"
          spacing={4}
          sx={{ width: '50%', alignItems: 'center' }}
        >

          <CircleIcon sx={{
            fontSize: 30, marginLeft: '3%', color,
          }}
          />
          <Typography>
            {name}
          </Typography>
          <Typography sx={{ color: typographyColor }}>
            {status}
          </Typography>
          <Typography sx={{ color: typographyColor }}>
            {team}
          </Typography>
        </Stack>
        <Box sx={{
          p: 0.1, bgcolor: 'background.paper', borderRadius: 3, width: 'flex', mr: 3, display: 'flex', m: 2, alignItems: 'center',
        }}
        >
          <Button
            variant="outlined"
            // size="small"
            sx={{
              borderRadius: 100,
              mr: 2,
            }}
            style={{ textTransform: 'none' }}
            onClick={handleOpen}
          >
            Add To Team
          </Button>
          <Modal
            open={addTeam}
            OnClose={handleClose}
          >
            <Box sx={{
              backgroundColor: 'background.paper',
              position: 'absolute',
              top: '30%',
              left: '50%',
              p: '5%',
              borderRadius: 5,
              transform: 'translate(-50%, -50%)',

            }}
            >
              <Typography variant="h6" sx={{ paddingBottom: 1 }}>
                Select A Team
              </Typography>
              <div>
                <FormControl variant="standard" sx={{ m: 1, minWidth: 120 }}>
                  <InputLabel id="demo-simple-select-standard-label">Team</InputLabel>
                  <Select
                    labelId="demo-simple-select-standard-label"
                    id="demo-simple-select-standard"
                    value={team}
                    onChange={handleChange}
                    label="Team"
                  >
                    <MenuItem value="">
                      <em>None</em>
                    </MenuItem>
                    <MenuItem value="Team 1">Team 1</MenuItem>
                    <MenuItem value="Team 2">Team 2</MenuItem>
                    <MenuItem value="Team 3">Team 3</MenuItem>
                  </Select>
                </FormControl>
              </div>
              <Button sx={{ paddingTop: 4 }} size="small" onClick={() => setAddTeam(false)}>Confirm</Button>
              <Button sx={{ paddingTop: 4 }} size="small" onClick={() => setAddTeam(false)}>Close</Button>
            </Box>
          </Modal>

          <Button
            variant="contained"
            color="secondary"
            sx={{
              borderRadius: 100,
              mr: 2,
            }}
            style={{ textTransform: 'none' }}
            onClick={() => history.push(`./genemapper/result/${projectId}`)}
            disabled={status !== 'DONE'}
          >
            See Results
          </Button>
          <Button>
            <DeleteOutlineIcon />
          </Button>
          {/* <CustomButton
            type="primary"
            sx={{
              outerHeight: '70%',
            }}
          >
            Create Mapping
          </CustomButton> */}

        </Box>
      </Stack>
    </Card>
  );
}
