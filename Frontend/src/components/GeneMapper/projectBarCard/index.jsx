import * as React from 'react';
import Card from '@mui/material/Card';
import Button from '@mui/material/Button';
import Typography from '@mui/material/Typography';
import { Box, Stack } from '@mui/material';
import CircleIcon from '@mui/icons-material/Circle';
import { useHistory } from 'react-router-dom';

export default function ProjectBarCard({ projectId, name, status }) {
  const history = useHistory();
  const [color, setColor] = React.useState(status === 'DONE' ? 'lightGreen' : status === 'IN PROGRESS' ? 'orange' : 'red');
  const [typographyColor, setTypographyColor] = React.useState(status === 'UPLOAD FAILED' ? 'red' : 'black');

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
          >
            Add To Team
          </Button>
          <Button
            variant="contained"
            color="secondary"
            sx={{
              borderRadius: 100,
            }}
            style={{ textTransform: 'none' }}
            onClick={() => history.push(`./genemapper/result/${projectId}`)}
          >
            See Results
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
